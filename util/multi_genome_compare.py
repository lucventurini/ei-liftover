#!/usr/bin/env python3

import argparse
import pyfaidx
import parasail
from Bio import Seq
import itertools
from Mikado.parsers.bed12 import Bed12Parser, BED12
from Mikado.utilities.overlap import overlap
import logging
import logging.handlers
import re
import numpy as np
from eiliftover.util.contrast import array_compare
from collections import defaultdict
import sys
import multiprocessing as mp
import sqlite3
import tempfile
import os
import operator
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.engine import create_engine
from sqlalchemy import Column, Integer, Float, Text, CHAR, Index
from sqlalchemy.orm.session import Session


DBBASE = declarative_base()

__doc__ = """"""


class BedIndex(DBBASE):

    __tablename__ = "bed"
    name = Column(Text, unique=True, primary_key=True)
    nidx = Index("sidx", "name", unique=True)
    start = Column(Integer)
    # end = Column(Integer)

    def __init__(self, name, start):

        self.name, self.start = name, start

        
def memoize_bed(string, sql):

    """Function to memoize a BED12 file for fast access"""

    db = sqlite3.connect(sql)
    db.execute("CREATE TABLE IF NOT EXISTS bed (chrom TEXT, start INT, \
    end INT, name TEXT, score REAL, strand CHAR, thick_start INT, \
    thick_end INT, rgb TEXT, count INT, sizes TEXT, starts TEXT)")
    db.execute("CREATE INDEX IF NOT EXISTS bedidx ON bed(name)")

    # records = dict()
    beds = []

    engine = create_engine("sqlite:///{}".format(sql))
    DBBASE.metadata.create_all(engine)

    session = Session(bind=engine, autocommit=True, autoflush=True)

    counter = 0
    with open(string, "rb") as parser:
        pos = parser.tell()  # This will be 0 as we are at the beginning
        for line in parser:
            fields = line.split(b"\t")
            if len(fields) != 12 or (fields and fields[0][0] == b"#"):
                header = True
            else:
                header = False

            if not header:
                name = fields[3].decode()
                if "ID=" in name:
                    groups = dict(re.findall("([^(;|=)]*)=([^;]*)", name))
                    name = groups["ID"]
                beds.append(BedIndex(name, pos))
                counter += 1
            pos += len(line)

            if len(beds) > 10 ** 6:
                session.begin(subtransactions=True)
                session.add_all(beds)
                session.commit()
                beds = []

    session.begin(subtransactions=True)
    session.add_all(beds)
    session.commit()
    session.close()
    res = engine.execute("select count(*) from bed").fetchone()
    assert res[0] == counter, (res[0], counter)


class ComparisonWorker(mp.Process):

    op_consumes = {
        "M": (True, True),
        "=": (True, True),
        "I": (True, False),
        "D": (False, True),
        "N": (False, True),
        "S": (True, False),
        "H": (False, False),
        "P": (False, False),
        "X": (True, True)
    }

    cigar_pattern = re.compile("({})".format("|".join(list(op_consumes.keys()))))

    def __init__(self, bed12records, bedfile, logging_queue, cdnas, entrance, identifier, consider_reference=False):

        super().__init__()
        self.logging_queue = logging_queue
        self.consider_reference = consider_reference
        self.handler = logging.handlers.QueueHandler(self.logging_queue)
        self.log = logging.getLogger(self.name)
        self.log.addHandler(self.handler)
        self.log.propagate = False
        self.log.setLevel("DEBUG")
        self.identifier = identifier
        self.name = "Comparer-{}".format(identifier)
        self.bed12records = sqlite3.connect(bed12records)
        self.__found_in_bed = dict()
        self.bedfile = open(bedfile, "rt")
        self.entrance = entrance
        self.cdnas = cdnas
        
        self.tmp_db_name = tempfile.NamedTemporaryFile(suffix="", prefix=".db", dir=os.getcwd())
        self.tmp_db = sqlite3.connect(self.tmp_db_name.name)
        self.tmp_db.execute("CREATE TABLE details (gid INT PRIMARY_KEY, row BLOB NOT NULL)")
        self.tmp_db.execute("CREATE TABLE summary (gid INT PRIMARY_KEY, row BLOB NOT NULL)")

    def prepare_info(self, transcript):

        cdna = str(self.fai[transcript]).upper()
        bed_position = self.__found_in_bed[transcript]
        self.bedfile.seek(bed_position)
        bed = BED12(self.bedfile.readline())
        assert bed.name == transcript, (bed.name, transcript, bed_position)
        bed = bed.to_transcriptomic(sequence=cdna)
        pep = str(Seq.Seq(str(cdna[max(0, bed.thick_start - 1):bed.thick_end])).translate())
        return cdna, bed, pep

    @classmethod
    def get_and_prepare_cigar(cls, t1cdna, t2cdna):

        # cigar_pattern = re.compile("(=|M|D|I|X|S|H)")
        result = parasail.sg_trace_scan_32(t1cdna, t2cdna, 11, 1, parasail.blosum62)
        # print(result.cigar.decode)
        values = re.split(cls.cigar_pattern, result.cigar.decode.decode())
        values = [(int(values[_ * 2]), values[_ * 2 + 1]) for _ in range(int((len(values) - 1) / 2))]
        return result, values

    @classmethod
    def create_translation_array(cls, cigar, query_exons, target_exons):

        """Inspired by https://github.com/MullinsLab/Bio-Cigar/blob/master/lib/Bio/Cigar.pm
        This function will translate the exons into a array-like space, allowing for
        comparison of the two structures.

        The main problem here is how we want to visualise


        """

        # Format: operation: [consumes_query, consumes_target]
        # First thing: determine the total length of the common space

        common_length = cls.cigar_length_in_common(cigar)
        # print(common_length)
        query_array = [None] * common_length
        target_array = [None] * common_length

        # Now we have to translate the exons into this common space
        # We will do this by creating "exons" derived from the alignment

        common_pos = 0
        query_pos = 0
        target_pos = 0

        for length, op in cigar:
            consumes_query, consumes_target = cls.op_consumes[op]
            if not any((consumes_query, consumes_target)):
                continue
            else:
                for pos in range(common_pos, common_pos + length):
                    if consumes_query:
                        query_pos += 1
                        query_array[pos] = query_pos
                    if consumes_target:
                        target_pos += 1
                        target_array[pos] = target_pos
                common_pos += length

        for pos in range(1, query_exons[-1][1] + 1):
            assert pos in query_array, (pos, min([_ for _ in query_array if _ is not None]), max([_ for _ in query_array if _ is not None]))
        for pos in range(1, target_exons[-1][1] + 1):
            assert pos in target_array

        try:
            c_query_exons = cls.transfer_exons(query_exons, query_array)
        except ValueError as exc:
            raise ValueError("\n".join([str(_) for _ in
                                        [query_exons, query_array, cigar, exc]]))
        try:
            c_target_exons = cls.transfer_exons(target_exons, target_array)
        except ValueError as exc:
            raise ValueError("\n".join([str(_) for _ in
                                        [target_exons, target_array, cigar, exc]]))

        return c_query_exons, c_target_exons, list(zip(query_array, target_array))

    @classmethod
    def cigar_length_in_common(cls, cigar):

        return sum(length for length, op in cigar if any(cls.op_consumes[op]))

    @staticmethod
    def transfer_exons(exons, c_array):
        """This static method will translate one set of exons into a different coordinate system
        (given an array that translates each position in system 1 to a position in system 2)"""

        c_exons = []
        for exon in exons:
            start, end = exon
            c_start = c_array.index(start)
            c_end = c_array.index(end)
            c_exons.append((c_start, c_end))

        return c_exons

    def run(self):

        self.__found_in_bed = dict(_ for _ in self.bed12records.execute("SELECT name, start from bed"))
        self.fai = pyfaidx.Fasta(self.cdnas)
        self.log.debug("Finished loading information for the process")
        
        while True:
            group, cases = self.entrance.get()
            if group == "EXIT":
                self.log.debug("EXIT for %s", self.name)
                # self.entrance.task_done()
                self.entrance.put(("EXIT", None))
                break
            else:
                details, summary = self.analyse_group(group, cases)
                if summary is not None:
                    for detail in details:
                        try:
                            row = "|".join(detail)
                        except TypeError:
                            raise TypeError(detail)
                        self.tmp_db.execute("INSERT INTO details VALUES (?, ?)", (group, row))

                    # [self.tmp_db.execute("INSERT INTO details VALUES (?, ?)", (group, "|".join(_))) for _ in details]
                    self.tmp_db.execute("INSERT INTO summary VALUES (?, ?)", (group, "|".join(summary)))
                    self.tmp_db.commit()
            continue
        self.tmp_db.commit()
        self.tmp_db.close()
        assert os.path.exists(self.tmp_db_name.name)

        return

    def _analyse_cDNAs(self, cdnas, beds, peps):

        result, cigar = self.get_and_prepare_cigar(*cdnas)

        t1bed, t2bed = beds

        assert sorted(t1bed.blocks)[-1][1] == len(cdnas[0]), (t1bed.blocks, len(cdnas[0]))
        assert sorted(t2bed.blocks)[-1][1] == len(cdnas[1]), (t2bed.blocks, len(cdnas[1]))

        try:
            c_t1_exons, c_t2_exons, common = self.create_translation_array(cigar, t1bed.blocks, t2bed.blocks)
        except (ValueError, AssertionError) as exc:
            raise ValueError(exc)

        # Common: list(zip(query_array, target_array))

        identical = sum(length for length, op in cigar if op in ("M", "="))
        if identical == 0:
            raise ValueError(cigar)
        identity = round(100 * identical / len(common), 2)
        result = array_compare(np.ravel(np.array(c_t1_exons)),
                               np.ravel(np.array(c_t2_exons)), identity)
        result, ccode = result[:-1].reshape((2, 3)), int(result[-1])
        # Now that we have analysed the cDNAs, it is time for the CDS

        if t1bed.coding and t2bed.coding:

            t1_coding_exons = [(max(t1bed.thick_start - 1, _[0]), min(t1bed.thick_end, _[1])) for _ in t1bed.blocks
                               if overlap(_, (t1bed.thick_start - 1, t1bed.thick_end)) > 0]
            assert t1_coding_exons, (
            t1bed.blocks, t1bed.block_starts, t1bed.block_sizes, t1bed.thick_start, t1bed.thick_end)
            t2_coding_exons = [(max(t2bed.thick_start - 1, _[0]), min(t2bed.thick_end, _[1])) for _ in t2bed.blocks
                               if overlap(_, (t2bed.thick_start - 1, t2bed.thick_end)) > 0]
            assert t2_coding_exons

            query_array, target_array = list(zip(*common))
            c_t1_coding = self.transfer_exons(t1_coding_exons, query_array)
            c_t2_coding = self.transfer_exons(t2_coding_exons, target_array)

            t1pep, t2pep = peps

            self.log.debug("CDS: %s:%s-%s: %s",
                           t2bed.chrom, t2bed.thick_start - 2, t2bed.thick_end, t2pep)
            # print(t2bed.chrom, t2bed.thick_start-2, t2bed.thick_end, t2_coding_exons, t2pep)

            coding_result, coding_cigar = self.get_and_prepare_cigar(t1pep, t2pep)
            coding_common = self.cigar_length_in_common(coding_cigar)
            coding_identical = sum(length for length, op in coding_cigar if op in ("M", "="))
            if coding_identical == 0:
                raise ValueError((coding_cigar, t1pep, t2pep))
            coding_identity = round(100 * coding_identical / coding_common, 2)

            coding_result = array_compare(np.ravel(np.array(c_t1_coding, dtype=np.int)),
                                          np.ravel(np.array(c_t2_coding, dtype=np.int)),
                                          coding_identity)
            coding_result, coding_ccode = coding_result[:-1].reshape((2, 3)), int(coding_result[-1])

        else:
            c_t1_coding = "NA"
            c_t2_coding = "NA"
            coding_identity = 0
            coding_result = np.zeros((2, 3))
            coding_ccode = 0

        return (c_t1_exons, c_t2_exons, identity, result, ccode,
                c_t1_coding, c_t2_coding, coding_identity, coding_result, coding_ccode)

    def _analyse_combination(self, t1, t2):

        cdnas, beds, peps = list(zip(t1, t2))  
        t1, t2 = beds[0].name, beds[1].name        
        
        # Get the results for the cDNAs
        # beds = [beds[0].to_transcriptomic(sequence=cdnas[0]), beds[1].to_transcriptomic(sequence=cdnas[1])]

        c_t1_exons, c_t2_exons, identity, result, ccode, \
            c_t1_coding, c_t2_coding, coding_identity, coding_result, coding_ccode = self._analyse_cDNAs(cdnas, beds, peps)
        # Now let's get the results for the proteins

        if ccode > 0:
            ccode = chr(ccode)
        else:
            ccode = "NA"
        if coding_ccode > 0:
            coding_ccode = chr(coding_ccode)
        else:
            coding_ccode = "NA"

        deta = (t1, str(c_t1_exons), t2, str(c_t2_exons), identity,  # c_t1_exons, c_t2_exons,
                *["{:0.2f}".format(100 * _) for _ in result[0]],
                *["{:0.2f}".format(100 * _) for _ in result[1]],
                ccode,
                c_t1_coding, c_t2_coding,
                coding_identity,
                *["{:0.2f}".format(100 * _) for _ in coding_result[0]],
                *["{:0.2f}".format(100 * _) for _ in coding_result[1]],
                coding_ccode)
        deta = (str(_) for _ in deta)

        return deta, result, identity, coding_result, coding_identity

    def analyse_group(self, group, cases):

        exon_f1 = []
        junc_f1 = []
        cds_f1 = []
        cds_junc_f1 = []
        iden = []
        cds_iden = []
        details, summary = [], None

        to_remove = set()
        if len(cases) < 2:
            self.log.error("Wrong case: %s", ",".join(cases))
            return details, summary

        self.log.debug("Group %s: cases %s", group, ", ".join(cases))
        for tid in cases:
            if tid not in self.fai or tid not in self.__found_in_bed:
                self.log.warning("%s not found among records for group %s, expunging.", tid, group)
                to_remove.add(tid)

        if len(cases) - len(to_remove) < 2:
            self.log.warning("Not enough records kept for group %s; removed: %s", group, ",".join(to_remove))
            return details, summary
        [cases.remove(_) for _ in to_remove]

        data = dict()
        for name in cases:
            data[name] = self.prepare_info(name)

        if self.consider_reference is True:
            combs = itertools.zip_longest([cases[0]], cases[1:], fillvalue=cases[0])
        else:
            combs = itertools.combinations(cases, 2)
        
        for comb in combs:
            self.log.debug("Combination: %s, %s", *comb)
            deta, result, identity, coding_result, coding_identity = self._analyse_combination(data[comb[0]], data[comb[1]])

            details.append(deta)
            exon_f1.append(result[0][2])
            junc_f1.append(result[1][2])
            cds_f1.append(coding_result[0][2])
            cds_junc_f1.append(coding_result[1][2])
            iden.append(identity)
            cds_iden.append(coding_identity)

        summary = ("{:0.2f}".format(min(100 * iden)),
                   "{:0.2f}".format(max(100 * iden)),
                   "{:0.2f}".format(min(junc_f1)),
                   "{:0.2f}".format(max(junc_f1)),
                   "{:0.2f}".format(min(exon_f1)),
                   "{:0.2f}".format(max(exon_f1)),
                   "{:0.2f}".format(min(100 * cds_iden)),
                   "{:0.2f}".format(max(100 * cds_iden)),
                   "{:0.2f}".format(min(cds_junc_f1)),
                   "{:0.2f}".format(max(cds_junc_f1)),
                   "{:0.2f}".format(min(cds_f1)),
                   "{:0.2f}".format(max(cds_f1)),
                   ",".join(cases))
        return details, summary


class OutPrinter(mp.Process):

    def __init__(self, name, sent, dbnames, logging_queue, summary=False):
        self.out_name, self.dbnames = name, dbnames
        self.summary = summary
        self.logging_queue = logging_queue
        super().__init__()
        self.sent = sent
        self.log = logging.getLogger(self._name)
        self.handler = logging.handlers.QueueHandler(self.logging_queue)
        self.log.addHandler(self.handler)
        self.log.propagate = False

    def print_summary(self):
        dbs = [sqlite3.connect(_.name) for _ in self.dbnames]
        rows = []
        with open(self.out_name, "wt") as out:
            header = "Group Min(Identity) Max(Identity)".split()
            header.extend(["Min(Junction F1)", "Max(Junction F1)"])
            header.extend(["Min(Exon F1)", "Max(Exon F1)"])
            header.extend("Min(CDS_Identity) Max(CDS_Identity)".split())
            header.extend(["Min(CDS Junction F1)", "Max(CDS Junction F1)"])
            header.extend(["Min(CDS Exon F1)", "Max(CDS Exon F1)"])
            header.append("IDs")
            print(*header, file=out, sep="\t")
            for db in dbs:
                rows.extend(db.execute("SELECT * FROM summary").fetchall())
                continue

            for row in sorted(rows, key=operator.itemgetter(0)):
                print(row[0], *row[1].split("|"), sep="\t", file=out)

        [_.close() for _ in dbs]
        return

    def print_detailed(self):
        dbs = [sqlite3.connect(_.name) for _ in self.dbnames]
        # rows = []
        header = "Group T1 T1_exons".split()
        header += "T2 T2_exons".split()
        header += ["Identity"]
        header += ["{}(Exon)".format(_) for _ in ["Recall", "Precision", "F1"]]
        header += ["{}(Junction)".format(_) for _ in ["Recall", "Precision", "F1"]]
        header += ["CCode"]
        header += ["T1_exons(CDS)", "T2_exons(CDS)"]
        header += ["Identity(CDS)"]
        header += ["{}(CDS,Exon)".format(_) for _ in ["Recall", "Precision", "F1"]]
        header += ["{}(CDS,Junction)".format(_) for _ in ["Recall", "Precision", "F1"]]
        header += ["CCode(CDS)"]

        with open(self.out_name, "wt") as detailed:
            print(
                *header,
                sep="\t", file=detailed)
            for index, group in enumerate(sorted(self.sent.keys())):
                db = dbs[self.sent[group]]
                details = db.execute("SELECT row FROM details WHERE gid=?", (str(group),)).fetchall()
                for detail in details:
                    print(group, *detail[0].split("|"), sep="\t", file=detailed)
                # if index in self.stopgaps:
                #     self.log.info("Done {}% of the summary", self.stopgaps.index(index) * 5)
        [_.close() for _ in dbs]
        return

    def run(self):
        if self.summary is False:
            self.print_detailed()
        else:
            self.print_summary()


def main():

    """"""

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--cdnas", required=True,
                        help="FASTA file with the cDNAs to be analysed.")
    parser.add_argument("-r", "--reference", action="store_true",
                        default=False,
                        help="""Flag. If set, the first element in each group is considered the 'reference'
                        and pairwise comparisons will be performed only against this model.""")
    parser.add_argument("--groups", required=True,
                        help="Tab separated file with the groups to be analysed.")
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument("-q", "--quiet", default=False, action="store_true")
    verbosity.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument("--bed12", required=True,
                        help="BED12 file with the coordinates of the models to analyse.")

    parser.add_argument("-d", "--detailed", type=argparse.FileType("wt"), required=True)
    parser.add_argument("-t", "--threads", type=int, default=mp.cpu_count())
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("wt"), required=True)
    args = parser.parse_args()

    # log = create_default_logger("log", level="INFO")

    # Step 1: memorise the BED12 for fast access
    # bed12records = memoize_bed(args.bed12)
    # Step 2: FAI of the cDNAs
    # Step 3: for each group in the groups file, perform a pairwise comparison

    groups = defaultdict(list)

    with open(args.groups) as group_file:
        for line in group_file:
            if line[0] == "#":
                continue
            tid, group = line.rstrip().split()
            groups[group].append(tid)

    formatter = logging.Formatter("{asctime} - {name} - {filename}:{lineno} - {levelname} - {funcName} \
- {processName} - {message}",
        style="{"
    )

    logging_queue = mp.Queue(-1)
    logger = logging.getLogger("listener")
    logger.propagate = False
    log_handler = logging.StreamHandler()
    if args.verbose:
        level = "DEBUG"
    elif args.quiet:
        level = "WARNING"
    else:
        level = "INFO"
    logger.setLevel(level)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)
    logger.info("Starting to load BED file")
    # bed_db = tempfile.NamedTemporaryFile(delete=True, suffix="bed_database", prefix=".db", dir=os.getcwd())
    bed_db = args.bed12 + ".bidx"
    if not os.path.exists(bed_db) or os.stat(bed_db).st_ctime < os.stat(args.bed12).st_ctime:
        # print(bed_db)
        logger.info("Loading ")
        memoize_bed(args.bed12, bed_db)
        logger.info("Loaded BED file")

    writer = logging.handlers.QueueListener(logging_queue, logger)
    writer.respect_handler_level = True
    writer.start()

    send_queue = mp.Queue(-1)

    procs = [ComparisonWorker(bed_db, args.bed12, logging_queue, args.cdnas, send_queue, _, consider_reference=args.reference)
             for _ in range(args.threads)]
    [_.start() for _ in procs]

    for group, cases in groups.items():
        send_queue.put((group, cases))

    send_queue.put(("EXIT", None))
    [_.join() for _ in procs]

    logger.info("Finished performing the direct comparisons, gathering info from the SQLite databases")
    dbnames = [_.tmp_db_name for _ in procs]
    logger.debug("DBs: {}".format(",".join([_.name for _ in dbnames])))
    dbs = [sqlite3.connect(_.name) for _ in dbnames]
    sent = dict()
    for pos, db in enumerate(dbs):
        for num in db.execute("SELECT gid FROM summary").fetchall():
            sent[num[0]] = pos
    [_.close() for _ in dbs]

    logger.info("Finished gathering info from the databases, starting to print")
    procs = [OutPrinter(name, sent, dbnames, logging_queue, summary=s) for name, s in zip((args.out.name, args.detailed.name),
                                                                                          (True, False))]
    [_.start() for _ in procs]
    [_.join() for _ in procs]
    logger.info("Finished")


main()
