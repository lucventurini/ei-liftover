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

    result = parasail.sg_trace_scan_16(str(t1cdna), str(t2cdna), 11, 1, parasail.blosum62)
    values = re.split(cigar_pattern, result.cigar.decode)
    values = [(int(values[_ * 2]), values[_ * 2 + 1]) for _ in range(int((len(values) - 1) / 2))]
    return result, values


def create_translation_arrat(cigar):

    """Inspired by https://github.com/MullinsLab/Bio-Cigar/blob/master/lib/Bio/Cigar.pm
    This function will create a double array, where each of the possible positions of query1 are associated with its
    positions on query 2.
    """

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

    q_2_t = dict()
    t_2_q = dict()

    start_target_pos, end_target_pos = 0, 0
    start_query_pos, end_query_pos = 0, 0

    for length, op in cigar:
        consumes_query, consumes_target = op_consumes[op]
        if consumes_query and consumes_target:

            end_query_pos += length
            end_target_pos += length
            # The two have the same length
            for query_key, target_key in itertools.zip_longest(range(start_query_pos, end_query_pos),
                                                    range(start_target_pos, end_target_pos)):
                q_2_t[query_key] = target_key
                t_2_q[target_key] = query_key
        elif consumes_target:
            # The query is skipped:
            end_target_pos += length
            for query_key, target_key in itertools.zip_longest([start_query_pos],
                                                               range(start_target_pos, end_target_pos),
                                                               fillvalue=start_query_pos):
                t_2_q[target_key] = query_key
            q_2_t[start_query_pos] = range(start_target_pos, end_target_pos)

        elif consumes_query:
            # The target is skipped:
            end_query_pos += length
            for query_key, target_key in itertools.zip_longest(range(start_query_pos, end_query_pos),
                                                               [start_target_pos],
                                                               fillvalue=start_target_pos):
                q_2_t[query_key] = target_key
            t_2_q[start_target_pos] = range(start_query_pos, end_query_pos)

        elif not consumes_query and not consumes_target:
            continue

    return q_2_t, t_2_q


def main():

    """"""

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--cdnas", required=True,
                        help="FASTA file with the cDNAs to be analysed.")
    parser.add_argument("--groups", required=True,
                        help="Tab separated file with the groups to be analysed.")
    parser.add_argument("--bed12", required=True,
                        help="BED12 file with the coordinates of the models to analyse.")
    args = parser.parse_args()

    log = create_default_logger("log", level="INFO")

    # Step 1: memorise the BED12 for fast access
    bed12records = memoize_bed(args.bed12)
    # Step 2: FAI of the cDNAs
    fai = pyfaidx.Fasta(args.cdnas)
    # Step 3: for each group in the groups file, perform a pairwise comparison

    with open(args.groups) as groups:
        for line in groups:
            cases = line.rstrip().split()
            for comb in itertools.combinations(cases, 2):
                t1, t2 = comb
                correct, cdnas, beds = prepare_info(t1, t2, fai, bed12records, log)
                if correct is False:
                    continue
                # Now the complicated part ...
                # TODO: Check proper alignment parameters
                t1cdna, t2cdna = cdnas
                result, cigar = get_and_prepare_cigar(t1cdna, t2cdna)

                q_2_t, t_2_q = create_translation_arrat(cigar)
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
                    # Here we are going to translate ...
                    exon = (t2_ar[pos], t2_ar[pos + 1] - 1)
                    trans_exon_start = t_2_q[exon[0]]
                    trans_exon_end = t_2_q[exon[1]]
                    if not isinstance(trans_exon_start, int):
                        # This is RANGE of values
                        trans_exon_start = int(trans_exon_start.start)

                    if not isinstance(trans_exon_end, int):
                        trans_exon_end = int(trans_exon_end.end)
                    t2_exons.append((int(trans_exon_start), int(trans_exon_end)))
                t1_transcript = Transcript()
                t1_transcript.chrom = "foo"
                t1_transcript.add_exons(t1_exons)
                t1_transcript.strand = "+"
                t1_transcript.id = t1

                t2_transcript = Transcript()
                t2_transcript.chrom = "foo"
                t2_transcript.add_exons(t2_exons)
                t2_transcript.strand = "+"
                t2_transcript.id = t2

                comp = Assigner.compare(t1_transcript, t2_transcript)
                print(comp)


    return


main()
