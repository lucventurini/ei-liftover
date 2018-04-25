#!/usr/bin/env python3


import sys
import argparse
from eiliftover.util import transfer_utilities as transfer
import re
import pyfaidx
from Mikado.parsers.bed12 import BED12, Bed12Parser
from Mikado.parsers.GFF import GFF3
from Mikado.transcripts import Transcript
from Bio import Seq
import parasail


__doc__ = """Script to try to translate the CDS from one coordinate system to another."""


def transfer_by_alignment(ref_pep, target_cdna, target_bed):
    frames = dict()
    # Get the three-frame translation
    for frame in range(3):
        frames[frame] = str(Seq.Seq(str(target_cdna[frame:])).translate(to_stop=False))

    # This will get the best match in the 3-frame translation
    aln, best_frame, score, best_cigar = None, None, float("-inf"), None
    for frame in frames:
        res, cigar = transfer.get_and_prepare_cigar(ref_pep, frames[frame], open=3, extend=1, matrix=parasail.blosum85)
        if res.score > score:
            aln, best_frame, score, best_cigar = res, frame, res.score, cigar

    # Now it is time to try to transfer it ... Ignore any deletions at the beginning
    cig_start = 0
    translation_start = 0

    while not transfer.op_consumes[best_cigar[cig_start][1]][0]:
        translation_start += best_cigar[cig_start][0]
        cig_start += 1

    # print(cig_start, best_cigar[cig_start:])
    # This is 0-based; we have to add 1 because we start 1 base after the gap at the beginning
    if translation_start > 0:
        translation_start = 3 * translation_start + 1

    translated = str(Seq.Seq(str(target_cdna[translation_start:])).translate(to_stop=True))
    # Logic to handle when the CDS is broken
    # This is 1-based, so we have to add 1 to
    target_bed.thick_start = translation_start + 1
    end = target_bed.thick_start + len(translated) * 3 - 1

    if translated and translated[0] != ref_pep[0]:
        if translation_start in (0, 1, 2):
            target_bed.phase = translation_start
        else:
            target_bed.coding = False
            return target_bed

    # Now check whether we can add the stop codon
    if end + 3 < len(target_cdna):
        end += 3
    else:  # Here we have to presume that it is open.
        end = len(target_cdna)

    # print(translation_start * 3, translated)
    target_bed.thick_end = end
    target_bed.coding = True
    return target_bed


def transfer_cds(transcript, ref_cdna, ref_bed, target_cdna, target_bed, logger=None, out=sys.stdout):

    if transcript is None:
        return

    transcript.finalize()

    orig_start, orig_end = target_bed.thick_start, target_bed.thick_end

    result, cigar = transfer.get_and_prepare_cigar(str(ref_cdna), str(target_cdna))
    ref_array, target_array = transfer.create_translation_array(cigar)

    ref_pep = str(Seq.Seq(str(ref_cdna[ref_bed.thick_start - 1:ref_bed.thick_end])).translate())

    # target_start, target_end = None, None

    try:
        target_start = target_array[ref_array.index(ref_bed.thick_start)]
    except IndexError:
        target_start = target_bed.start
    try:
        target_end = target_array[ref_array.index(ref_bed.thick_end)]
    except IndexError:
        target_end = target_bed.end

    # Let us check that the protein is functional ...

    replace = True

    if target_start is not None and target_end is not None:
        if target_start == target_bed.thick_start and target_end == target_bed.thick_end:
            replace = False  # Everything fine here
        else:
            target_pep = str(Seq.Seq(str(target_cdna[target_start - 1:target_end])).translate())
            if (target_pep[0] == ref_pep[0] and target_pep[-1] == ref_pep[-1] and
                    ((ref_pep[-1] == "*" and target_pep.count("*") == 1) or target_pep.count("*") == 0)):
                target_bed.thick_start, target_bed.thick_end = target_start, target_end
            else:
                target_bed.coding = False
    else:
        target_bed.coding = False

    if not replace:
        transcript.attributes["original_cds"] = True
        transcript.attributes["aligner_cds"] = True
        print("CDS transferred correctly for {}".format(ref_bed.chrom))
    else:
        transcript.strip_cds()
        transcript.attributes["original_cds"] = True
        transcript.attributes["aligner_cds"] = False
        target_bed.transcriptomic = True
        if target_bed.coding is True:
            try:
                valid = (not target_bed.invalid)
            except:
                valid = False
        else:
            valid = False
        if valid:
            print("Transferred original CDS for {}: from {}-{} to {}-{}".format(ref_bed.chrom,
                                                               orig_start, orig_end,
                                                               target_bed.thick_start, target_bed.thick_end))

            transcript.load_orfs([target_bed])

        else:
            # Here we try to transfer the CDS through alignment
            transcript.attributes["original_cds"] = False
            transcript.attributes["original_cds_coords"] = "{}-{}".format(target_start, target_end)
            recalc_coords = (target_bed.thick_start, target_bed.thick_end)
            target_bed = transfer_by_alignment(ref_pep, target_cdna, target_bed)
            target_bed.transcriptomic = True
            # Let's see what happens when we just use the transferred cDNA start

            if target_bed.coding is True and target_bed.invalid is False:
                if (orig_start, orig_end) == (target_bed.thick_start, target_bed.thick_end):
                    print("GMAP corrected the CDS for {}: from {}-{} to {}-{}".format(
                        ref_bed.chrom,
                        recalc_coords[0], recalc_coords[1],
                        target_bed.thick_start, target_bed.thick_end
                    ))
                    transcript.attributes["aligner_cds"] = True
                else:
                    print("New CDS for {}: from {}-{} to {}-{}".format(ref_bed.chrom,
                                                                       orig_start, orig_end,
                                                                       target_bed.thick_start, target_bed.thick_end))
                transcript.load_orfs([target_bed])
            else:
                print("Stripping CDS from {}".format(ref_bed.chrom))
                # transcript.load_orfs([target_bed])

    print(transcript.format("gff3"), file=out)


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--bed12", nargs=2, required=True, help="Transcriptomic cDNAs BED12s")
    parser.add_argument("--cdnas", nargs=2, required=True)
    parser.add_argument("-gf", help="GFF3 of the transferred annotation.", required=True)
    parser.add_argument("--out", default=sys.stdout, type=argparse.FileType("wt"))
    args = parser.parse_args()

    cdnas = dict()
    beds = dict()

    gmap_pat = re.compile("\.mrna[0-9]*$")

    for key, cdna, bed in zip(("ref", "target"), args.cdnas, args.bed12):
        cdnas[key] = pyfaidx.Fasta(cdna)
        beds[key] = dict()
        for entry in Bed12Parser(bed):
            if entry.header:
                continue
            beds[key][re.sub(gmap_pat, "", entry.chrom)] = entry

    # Now let us start parsing the GFF3, which we presume being a GMAP GFF3
    transcript = None

    for line in GFF3(args.gf):
        if line.header is True or line.gene is True:
            print(line, file=args.out)
            continue
        elif line.is_transcript is True:
            if transcript:
                tid = re.sub(gmap_pat, "", transcript.id)
                transfer_cds(transcript,
                             cdnas["ref"][tid], beds["ref"][tid],
                             cdnas["target"][transcript.id], beds["target"][tid],
                             out=args.out)
            transcript = Transcript(line)
        elif line.is_exon is True:
            transcript.add_exon(line)

    if transcript:
        tid = re.sub(gmap_pat, "", transcript.id)
        transfer_cds(transcript,
                     cdnas["ref"][tid], beds["ref"][tid],
                     cdnas["target"][transcript.id], beds["target"][tid],
                     out=args.out)

    return


if __name__ == "__main__":
    main()
