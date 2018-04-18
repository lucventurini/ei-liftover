#!/usr/bin/env python3

import argparse
import pyfaidx
import parasail
import edlib
from Bio import Seq
import itertools
from Mikado.parsers.bed12 import Bed12Parser, BED12
from Mikado.utilities.overlap import overlap
from Mikado.utilities.log_utils import create_default_logger
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


__doc__ = """"""


def memoize_bed(string, sql):

    """Function to memoize a BED12 file for fast access"""

    db = sqlite3.connect(sql)
    db.execute("CREATE TABLE IF NOT EXISTS bed (chrom TEXT, start INT, \
    end INT, name TEXT, score REAL, strand CHAR, thick_start INT, \
    thick_end INT, rgb TEXT, count INT, sizes TEXT, starts TEXT)")
    db.execute("CREATE INDEX IF NOT EXISTS bedidx ON bed(name)")

    # records = dict()
    for record in Bed12Parser(string):
        if record.header is True:
            continue
        try:
            db.execute("INSERT INTO bed VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                       (
                           record.chrom, record.start, record.end, record.name, record.score,
                           record.strand, record.thick_start, record.thick_end, record.rgb,
                           record.block_count, ",".join(str(_) for _ in record.block_sizes),
                           ",".join(str(_) for _ in record.block_starts)
                        ))
        except sqlite3.InterfaceError as exc:
            raise sqlite3.InterfaceError("{}:\n{}".format(exc,
                                                          (
                                                              record.chrom, record.start, record.end, record.name,
                                                              record.score,
                                                              record.strand, record.thick_start, record.thick_end,
                                                              record.rgb,
                                                              record.block_count, record.block_sizes,
                                                              record.block_starts
                                                          )
                                                          ))
        db.commit()
        # records[record.name] = record
    db.close()
    # return records


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

    def __init__(self, bed12records, logging_queue, cdnas, entrance, identifier, consider_reference=False):

        super().__init__()
        self.logging_queue = logging_queue
        self.consider_reference = consider_reference
        self.handler = logging.handlers.QueueHandler(self.logging_queue)
        self.log = logging.getLogger(self.name)
        self.log.addHandler(self.handler)
        self.log.propagate = False
        self.identifier = identifier
        self.name = "Comparer-{}".format(identifier)
        self.bed12records = sqlite3.connect(bed12records)
        self.__found_in_bed = set([_[0] for _ in self.bed12records.execute("SELECT name from bed")])
        self.fai = pyfaidx.Fasta(cdnas)
        self.entrance = entrance

        self.tmp_db_name = tempfile.mktemp("", ".db", os.getcwd())
        if os.path.exists(self.tmp_db_name):
            os.remove(self.tmp_db_name)
        self.tmp_db = sqlite3.connect(self.tmp_db_name)
        self.tmp_db.execute("CREATE TABLE details (gid INT PRIMARY_KEY, row BLOB NOT NULL)")
        self.tmp_db.execute("CREATE TABLE summary (gid INT PRIMARY_KEY, row BLOB NOT NULL)")

    def prepare_info(self, t1, t2):

        t1cdna = str(self.fai[t1]).upper()
        t2cdna = str(self.fai[t2]).upper()
        t1bed = BED12("\t".join([str(_) for _ in
                                 self.bed12records.execute("SELECT * FROM bed WHERE name=?", (t1, )).fetchone()]))


        t2bed = BED12("\t".join([str(_) for _ in
                                 self.bed12records.execute("SELECT * FROM bed WHERE name=?", (t2, )).fetchone()]))

        return (t1cdna, t2cdna), (t1bed, t2bed)

    @classmethod
    def get_and_prepare_cigar(cls, t1cdna, t2cdna):

        # cigar_pattern = re.compile("(=|M|D|I|X|S|H)")
        result = parasail.sg_trace_scan_32(t1cdna, t2cdna, 11, 1, parasail.blosum100)
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
        query_pos = -1
        target_pos = -1

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

        try:
            c_query_exons = cls.transfer_exons(query_exons, query_array)
        except ValueError as exc:
            raise ValueError("\n".join([str(_) for _ in
                                        [query_exons, query_array]]))
        try:
            c_target_exons = cls.transfer_exons(target_exons, target_array)
        except ValueError as exc:
            raise ValueError("\n".join([str(_) for _ in
                                        [target_exons, target_array]]))

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
            c_start, c_end = c_array.index(start), c_array.index(end)
            c_exons.append((c_start, c_end))

        return c_exons

    def run(self):

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
        assert os.path.exists(self.tmp_db_name)

        return

    def _analyse_cDNAs(self, cdnas, beds):

        result, cigar = self.get_and_prepare_cigar(*cdnas)
        t1bed, t2bed = beds

        print(t1bed.blocks, t2bed.blocks)

        try:
            c_t1_exons, c_t2_exons, common = self.create_translation_array(cigar, t1bed.blocks, t2bed.blocks)
        except ValueError as exc:
            print(t1bed.blocks)
            print(t2bed.blocks)
            print(cigar)
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

        t1_coding_exons = [(max(t1bed.thick_start, _[0]), min(t1bed.thick_end, _[1])) for _ in t1bed.blocks
                           if overlap(_, (t1bed.thick_start, t1bed.thick_end)) > 0]
        assert t1_coding_exons, (t1bed.blocks, t1bed.block_starts, t1bed.block_sizes, t1bed.thick_start, t1bed.thick_end)
        t2_coding_exons = [(max(t2bed.thick_start, _[0]), min(t2bed.thick_end, _[1])) for _ in t2bed.blocks
                           if overlap(_, (t1bed.thick_start, t1bed.thick_end)) > 0]
        assert t2_coding_exons

        query_array, target_array = list(zip(*common))
        c_t1_coding = self.transfer_exons(t1_coding_exons, query_array)
        c_t2_coding = self.transfer_exons(t2_coding_exons, target_array)

        t1pep = str(Seq.Seq(str(cdnas[0][t1bed.thick_start-1:t1bed.thick_end + 1])).translate())
        t2pep = str(Seq.Seq(str(cdnas[1][t2bed.thick_start-1:t2bed.thick_end + 1])).translate())

        coding_result, coding_cigar = self.get_and_prepare_cigar(t1pep, t2pep)
        coding_common = self.cigar_length_in_common(coding_cigar)
        coding_identical = sum(length for length, op in coding_cigar if op in ("M", "="))
        if coding_identical == 0:
            raise ValueError((coding_cigar, t1pep, t2pep))
        coding_identity = round(100 * coding_identical / coding_common, 2)

        print(t1_coding_exons, t2_coding_exons, c_t1_coding, c_t2_coding, sep="\n")

        coding_result = array_compare(np.ravel(np.array(c_t1_coding, dtype=np.int)),
                                      np.ravel(np.array(c_t2_coding, dtype=np.int)),
                                      coding_identity)
        coding_result, coding_ccode = coding_result[:-1].reshape((2, 3)), int(result[-1])

        return (c_t1_exons, c_t2_exons, identity, result, ccode, coding_identity, coding_result, coding_ccode)

    def _analyse_combination(self, t1, t2):

        cdnas, beds = self.prepare_info(t1, t2)
        # Get the results for the cDNAs
        beds = [beds[0].to_transcriptomic(sequence=cdnas[0]), beds[1].to_transcriptomic(sequence=cdnas[1])]

        c_t1_exons, c_t2_exons, identity, result, ccode, \
            coding_identity, coding_result, coding_ccode = self._analyse_cDNAs(cdnas, beds)
        # Now let's get the results for the proteins

        deta = (t1, str(c_t1_exons), t2, str(c_t2_exons), identity, coding_identity, # c_t1_exons, c_t2_exons,
                *["{:0.2f}".format(100 * _) for _ in result[0]],
                *["{:0.2f}".format(100 * _) for _ in result[1]],
                chr(ccode), chr(coding_ccode))
        deta = (str(_) for _ in deta)

        return deta, result, identity

    def analyse_group(self, group, cases):

        exon_f1 = []
        junc_f1 = []
        iden = []
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

        if self.consider_reference is True:
            combs = itertools.zip_longest([cases[0]], cases[1:], fillvalue=cases[0])
        else:
            combs = itertools.combinations(cases, 2)

        for comb in combs:
            self.log.debug("Combination: %s, %s", *comb)
            deta, result, identity = self._analyse_combination(*comb)

            details.append(deta)
            exon_f1.append(result[0][2])
            junc_f1.append(result[1][2])
            iden.append(identity)

        summary = ("{:0.2f}".format(min(100 * iden)),
                   "{:0.2f}".format(max(iden)),
                   "{:0.2f}".format(min(junc_f1)),
                   "{:0.2f}".format(max(junc_f1)),
                   "{:0.2f}".format(min(exon_f1)),
                   "{:0.2f}".format(max(exon_f1)),
                   ",".join(cases))
        return details, summary


class OutPrinter(mp.Process):

    def __init__(self, name, dbnames, logging_queue, summary=False):
        self.out_name, self.dbnames = name, dbnames
        self.summary = summary
        self.logging_queue = logging_queue
        self.handler = logging.handlers.QueueHandler(self.logging_queue)
        self.log = logging.getLogger(self.name)
        self.log.addHandler(self.handler)
        self.log.propagate = False
        super().__init__()

    def print_summary(self):
        dbs = [sqlite3.connect(_) for _ in self.dbnames]
        rows = []
        with open(self.out_name, "wt") as out:
            header = "Group Min(Identity) Max(Identity)".split()
            header.extend(["Min(Junction F1)", "Max(Junction F1)"])
            header.extend(["Min(Exon F1)", "Max(Exon F1)"])
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
        dbs = [sqlite3.connect(_) for _ in self.dbnames]
        rows = []
        with open(self.out_name, "wt") as detailed:
            print(
                *"Group T1 T1_exons T2 T2_exons Identity Identity(CDS) Recall(Exon) Precision(Exon) F1(Exon) Recall(Junction) Precision(Junction) F1(Junction) CCode CCode(CDS)".split(),
                sep="\t", file=detailed)
            for index, group in enumerate(sorted(self.sent.keys())):
                db = dbs[self.sent[group]]
                details = db.execute("SELECT row FROM details WHERE gid=?", (str(group),)).fetchall()
                for detail in details:
                    print(group, *detail[0].split("|"), sep="\t", file=detailed)
                if index in self.stopgaps:
                    self.log.info("Done {}% of the summary", self.stopgaps.index(index) * 5)
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
    parser.add_argument("--bed12", required=True,
                        help="BED12 file with the coordinates of the models to analyse.")
    parser.add_argument("-d", "--detailed", type=argparse.FileType("wt"), required=True)
    parser.add_argument("-t", "--threads", type=int, default=mp.cpu_count())
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("wt"), required=True)
    args = parser.parse_args()

    log = create_default_logger("log", level="INFO")

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
    logger.setLevel("INFO")
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)
    logger.info("Starting to load BED file")
    bed_db = tempfile.mktemp(suffix="bed_database", prefix=".db", dir=os.getcwd())
    memoize_bed(args.bed12, bed_db)
    logger.info("Loaded BED file")

    writer = logging.handlers.QueueListener(logging_queue, logger)
    writer.start()

    send_queue = mp.Queue(-1)

    procs = [ComparisonWorker(bed_db, logging_queue, args.cdnas, send_queue, _, consider_reference=args.reference)
             for _ in range(args.threads)]
    [_.start() for _ in procs]

    for group, cases in groups.items():
        send_queue.put((group, cases))

    send_queue.put(("EXIT", None))
    [_.join() for _ in procs]

    logger.info("Finished performing the direct comparisons, gathering info from the SQLite databases")
    dbnames = [_.tmp_db_name for _ in procs]
    logger.debug("DBs: {}".format(",".join(dbnames)))
    dbs = [sqlite3.connect(_) for _ in dbnames]
    sent = dict()
    for pos, db in enumerate(dbs):
        for num in db.execute("SELECT gid FROM summary").fetchall():
            sent[num[0]] = pos
    [_.close() for _ in dbs]

    logger.info("Finished gathering info from the databases, starting to print")
    procs = [OutPrinter(name, sent, dbnames, logging_queue, summary=s) for name, s in zip((args.out, args.detailed),
                                                                                          (True, False))]
    [_.start() for _ in procs]
    [_.join() for _ in procs]
    logger.info("Finished")


main()
