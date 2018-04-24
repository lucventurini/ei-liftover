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


__doc__ = """Script to try to translate the CDS from one coordinate system to another."""


def transfer_cds(transcript, ref_cdna, ref_bed, target_cdna, target_bed, logger=None, out=sys.stdout):

    if transcript is None:
        return

    transcript.finalize()

    orig_start, orig_end = target_bed.thick_start, target_bed.thick_end

    result, cigar = transfer.get_and_prepare_cigar(str(ref_cdna), str(target_cdna))
    ref_array, target_array = transfer.create_translation_array(cigar)

    # ref_pep = ref_cdna[ref_bed.thick_start - 1:ref_bed.thick_end]

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
            target_pep = Seq.Seq(str(target_cdna[target_start - 1:target_end])).translate()
            if str(target_pep).startswith("M") and str(target_pep).endswith("*") and str(target_pep).count("*") == 1:
                target_bed.thick_start, target_bed.thick_end = target_start, target_end
            else:
                target_bed.coding = False
    else:
        target_bed.coding = False

    if not replace:
        print("CDS transferred correctly for {}".format(ref_bed.chrom))
    else:
        transcript.strip_cds()
        valid = True
        if target_bed.coding is True:
            try:
                valid = (not target_bed.invalid)
            except:
                valid = False
        else:
            valid = False
        if valid:
            print("New CDS for {}: from {}-{} to {}-{}".format(ref_bed.chrom,
                                                               orig_start, orig_end,
                                                               target_bed.thick_start, target_bed.thick_end))
            transcript.load_orfs([target_bed])
        else:
            print("Stripped CDS from {}".format(ref_bed.chrom))

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
