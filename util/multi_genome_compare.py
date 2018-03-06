#!/usr/bin/env python3

import argparse
import pyfaidx
import parasail
import edlib
import itertools
from Mikado.parsers.bed12 import Bed12Parser, BED12
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


__doc__ = """"""


def memoize_bed(string):

    """Function to memoize a BED12 file for fast access"""

    records = dict()
    for record in Bed12Parser(string):
        if record.header is True:
            continue
        records[record.name] = record

    return records


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
        self.bed12records = bed12records
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
        t1bed = self.bed12records[t1]
        t2bed = self.bed12records[t2]
        assert isinstance(t1bed, BED12)
        assert isinstance(t2bed, BED12)

        return (t1cdna, t2cdna), (t1bed, t2bed)

    @classmethod
    def get_and_prepare_cigar(cls, t1cdna, t2cdna):

        # cigar_pattern = re.compile("(=|M|D|I|X|S|H)")
        result = parasail.sg_trace_scan_32(t1cdna, t2cdna, 11, 1, parasail.blosum100)
        # print(result.cigar.decode)
        values = re.split(cls.cigar_pattern, result.cigar.decode)
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
        common_length = sum(length for length, op in cigar if any(cls.op_consumes[op]))

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

        c_query_exons = []
        c_target_exons = []

        for exon in query_exons:
            # Series of tuples
            start, end = exon
            c_start, c_end = query_array.index(start), query_array.index(end)
            c_query_exons.append((c_start, c_end))

        for exon in target_exons:
            start, end = exon
            c_start, c_end = target_array.index(start), target_array.index(end)
            c_target_exons.append((c_start, c_end))

        return c_query_exons, c_target_exons, list(zip(query_array, target_array))

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

    def _analyse_combination(self, t1, t2):

        cdnas, beds = self.prepare_info(t1, t2)
        result, cigar = self.get_and_prepare_cigar(*cdnas)
        self.log.debug("%s vs %s, CIGAR: %s", t1, t2, cigar)
        t1bed, t2bed = beds
        if t1bed.strand == "-":
            t1blocks = [0] + list(reversed(t1bed.block_sizes))
        else:
            t1blocks = [0] + t1bed.block_sizes
        if t2bed.strand == "-":
            t2blocks = [0] + list(reversed(t2bed.block_sizes))
        else:
            t2blocks = [0] + t2bed.block_sizes
        t1_ar = np.cumsum(t1blocks)
        t2_ar = np.cumsum(t2blocks)
        t1_exons = []
        t2_exons = []
        for pos in range(t1bed.block_count):
            t1_exons.append((int(t1_ar[pos]), int(t1_ar[pos + 1] - 1)))
        for pos in range(t2bed.block_count):
            t2_exons.append((int(t2_ar[pos]), int(t2_ar[pos + 1] - 1)))

        c_t1_exons, c_t2_exons, common = self.create_translation_array(cigar, t1_exons, t2_exons)

        identical = sum(length for length, op in cigar if op in ("M", "="))
        if identical == 0:
            raise ValueError(cigar)
        identity = round(100 * identical / len(common), 2)
        result = array_compare(np.ravel(np.array(c_t1_exons)),
                               np.ravel(np.array(c_t2_exons)), identity)
        result, ccode = result[:-1].reshape((2, 3)), int(result[-1])
        deta = (t1, str(c_t1_exons), t2, str(c_t2_exons), identity,  # c_t1_exons, c_t2_exons,
                *["{:0.2f}".format(100 * _) for _ in result[0]],
                *["{:0.2f}".format(100 * _) for _ in result[1]],
                chr(ccode))
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
            if tid not in self.fai or tid not in self.bed12records:
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
    bed12records = memoize_bed(args.bed12)
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
    logger.setLevel("WARNING")
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

    writer = logging.handlers.QueueListener(logging_queue, logger)
    writer.start()

    send_queue = mp.Queue(-1)

    procs = [ComparisonWorker(bed12records, logging_queue, args.cdnas, send_queue, _, consider_reference=args.reference)
             for _ in range(args.threads)]
    [_.start() for _ in procs]

    for group, cases in groups.items():
        send_queue.put((group, cases))

    send_queue.put(("EXIT", None))
    [_.join() for _ in procs]

    dbs = [_.tmp_db_name for _ in procs]
    logger.debug("DBs: {}".format(",".join(dbs)))
    dbs = [sqlite3.connect(_) for _ in dbs]
    sent = dict()
    for pos, db in enumerate(dbs):
        for num in db.execute("SELECT gid FROM summary").fetchall():
            sent[num[0]] = pos
    print(sent)

    print(
        *"Group T1 T1_exons T2 T2_exons Identity Recall(Exon) Precision(Exon) F1(Exon) Recall(Junction) Precision(Junction) F1(Junction) CCode".split(),
        sep="\t", file=args.detailed)
    header = "Group Min(Identity) Max(Identity)".split()
    header.extend(["Min(Junction F1)", "Max(Junction F1)"])
    header.extend(["Min(Exon F1)", "Max(Exon F1)"])
    header.append("IDs")
    print(*header, file=args.out, sep="\t")

    for group in groups:
        group = int(group)
        if group not in sent:
            log.warning("Group %s has been lost", group)
            continue
        db = dbs[sent[group]]
        details = db.execute("SELECT row FROM details WHERE gid=?", str(group)).fetchall()
        for detail in details:
            print(group, *detail[0].split("|"), sep="\t", file=args.detailed)
        summary = db.execute("SELECT row FROM details WHERE gid=?", str(group)).fetchone()
        print(group, *summary[0].split("|"), sep="\t", file=args.out)


main()
