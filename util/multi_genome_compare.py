#!/usr/bin/env python3

import argparse
import pyfaidx
import parasail
import edlib
import itertools
from Mikado.parsers.bed12 import Bed12Parser, BED12
from Mikado.utilities.log_utils import create_default_logger
import re
import numpy as np
from Mikado.scales.assigner import Assigner 
from Mikado.transcripts import Transcript
from eiliftover.util.contrast import array_compare
from collections import defaultdict
import sys


__doc__ = """"""


cigar_pattern = re.compile("(=|M|D|I|X|S|H)")


def memoize_bed(string):

    """Function to memoize a BED12 file for fast access"""

    records = dict()
    for record in Bed12Parser(string):
        if record.header is True:
            continue
        records[record.name] = record

    return records


def prepare_info(t1, t2, fai, bed12records, log):

    if t1 not in fai:
        log.warn("%s not in cDNA file, continuing", t1)
        return False, (), ()
    elif t2 not in fai:
        log.warn("%s not in cDNA file, continuing", t2)
        return False, (), ()
    t1cdna = fai[t1]
    t2cdna = fai[t2]
    if t1 not in bed12records:
        log.warn("%s not in BED12, continuing", t1)
        return False, (), ()
    elif t2 not in bed12records:
        log.warn("%s not in BED12, continuing", t2)
        return False, (), ()
    t1bed = bed12records[t1]
    t2bed = bed12records[t2]
    assert isinstance(t1bed, BED12)
    assert isinstance(t2bed, BED12)

    return True, (t1cdna, t2cdna), (t1bed, t2bed)


def get_and_prepare_cigar(t1cdna, t2cdna):

    # cigar_pattern = re.compile("(=|M|D|I|X|S|H)")
    result = parasail.sg_trace_scan_32(str(t1cdna).upper(), str(t2cdna).upper(), 11, 1, parasail.blosum100)
    # print(result.cigar.decode)
    values = re.split(cigar_pattern, result.cigar.decode)
    values = [(int(values[_ * 2]), values[_ * 2 + 1]) for _ in range(int((len(values) - 1) / 2))]
    return result, values


def create_translation_array(cigar, query_exons, target_exons):

    """Inspired by https://github.com/MullinsLab/Bio-Cigar/blob/master/lib/Bio/Cigar.pm
    This function will translate the exons into a array-like space, allowing for
    comparison of the two structures.

    The main problem here is how we want to visualise


    """

    # Format: operation: [consumes_query, consumes_target]
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

    # First thing: determine the total length of the common space
    common_length = sum(length for length, op in cigar if any(op_consumes[op]))

    # print(common_length)
    query_array = [None] * common_length
    target_array = [None] * common_length

    # Now we have to translate the exons into this common space
    # We will do this by creating "exons" derived from the alignment

    common_pos = 0
    query_pos = -1
    target_pos = -1

    for length, op in cigar:
        consumes_query, consumes_target = op_consumes[op]
        if not any((consumes_query, consumes_target)):
            continue
        else:
            for pos in range(common_pos, common_pos + length):
                if consumes_query:
                    query_pos +=1
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


def main():

    """"""

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--cdnas", required=True,
                        help="FASTA file with the cDNAs to be analysed.")
    parser.add_argument("--groups", required=True,
                        help="Tab separated file with the groups to be analysed.")
    parser.add_argument("--bed12", required=True,
                        help="BED12 file with the coordinates of the models to analyse.")
    parser.add_argument("-d", "--detailed", type=argparse.FileType("wt"), required=True)
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("wt"), required=True)
    args = parser.parse_args()

    log = create_default_logger("log", level="INFO")

    # Step 1: memorise the BED12 for fast access
    bed12records = memoize_bed(args.bed12)
    # Step 2: FAI of the cDNAs
    fai = pyfaidx.Fasta(args.cdnas)
    # Step 3: for each group in the groups file, perform a pairwise comparison

    groups = defaultdict(list)

    with open(args.groups) as group_file:
        for line in group_file:
            if line[0] == "#":
                continue
            tid, group = line.rstrip().split()
            groups[group].append(tid)

    print(*"Group T1 T1 exons T2 T2 exons Identity Recall(Exon) Precision(Exon) F1(Exon) Recall(Junction) Precision(Junction) F1(Junction) CCode".split(),
          sep="\t", file=args.detailed)
    header = "Group Min(Identity) Max(Identity)".split()
    header.extend(["Min(Junction F1)", "Max(Junction F1)"])
    header.extend(["Min(Exon F1)", "Max(Exon F1)"])
    header.append("IDs")
    print(*header, file=args.out, sep="\t")

    for group, cases in groups.items():

        exon_f1 = []
        junc_f1 = []
        iden = []

        for comb in itertools.combinations(cases, 2):
            t1, t2 = comb
            correct, cdnas, beds = prepare_info(t1, t2, fai, bed12records, log)
            if correct is False:
                continue
            # Now the complicated part ...
            # TODO: Check proper alignment parameters
            t1cdna, t2cdna = cdnas
            result, cigar = get_and_prepare_cigar(t1cdna, t2cdna)

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
            c_t1_exons, c_t2_exons, common = create_translation_array(cigar, t1_exons, t2_exons)

            identical = sum(length for length, op in cigar if op in ("M", "="))
            if identical == 0:
                raise ValueError(cigar)
            identity = round(100 * identical / len(common), 2)

            result = array_compare(np.ravel(np.array(c_t1_exons)),
                                   np.ravel(np.array(c_t2_exons)), identity)

            result, ccode = result[:-1].reshape((2, 3)), result[-1]

            print(group, t1, c_t1_exons, t2, c_t2_exons, identity,  # c_t1_exons, c_t2_exons,
                  *["{:0.2f}".format(100 * _) for _ in result[0]],
                  *["{:0.2f}".format(100 * _) for _ in result[1]],
                  ccode,
                  sep="\t", file=args.detailed)
            exon_f1.append(result[0][2])
            junc_f1.append(result[1][2])
            iden.append(identity)

        print(group,
              "{:0.2f}".format(min(100 * iden)),
              "{:0.2f}".format(max(iden)),
              "{:0.2f}".format(min(junc_f1)),
              "{:0.2f}".format(max(junc_f1)),
              "{:0.2f}".format(min(exon_f1)),
              "{:0.2f}".format(max(exon_f1)),
              ",".join(cases), sep="\t", file=args.out)


    return


main()
