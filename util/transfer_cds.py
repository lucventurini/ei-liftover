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
import multiprocessing as mp
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Integer, Column, BLOB
from sqlalchemy.orm import sessionmaker
import tempfile
import os


__doc__ = """Script to try to translate the CDS from one coordinate system to another."""

transfer_base = declarative_base()


class _Storer(transfer_base):

    __tablename__ = "storer"

    id = Column(Integer, primary_key=True)
    bed = Column(BLOB)
    gff3 = Column(BLOB)

    def __init__(self, id, bed, gff3):

        self.id, self.bed, self.gff3 = id, bed, gff3



class Transferer(mp.Process):

    def __init__(self, out_sq, queue):

        super().__init__()
        self.out_sq = out_sq
        self.engine = sqlalchemy.create_engine("sqlite:///{}".format(out_sq))
        transfer_base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        self.queue = queue

    def run(self):

        while True:
            obj_input = self.queue.get()
            if obj_input == "EXIT":
                self.queue.put("EXIT")
                self.session.commit()
                # self.join()
                return
            num_id, (transcript,
             ref_cdna, ref_bed,
             target_cdna, target_bed) = obj_input

            transcript, target_bed, pep_coords = transfer_cds(transcript, ref_cdna, ref_bed, target_cdna, target_bed)
            target_bed.name = "ID={};coding={}".format(target_bed.name, target_bed.coding)
            if target_bed.coding:
                target_bed.name = target_bed.name + ";phase={};original_pep_coords={}-{};original_pep_complete={}".format(
                    target_bed.phase, *pep_coords)
            else:
                target_bed.thick_start, target_bed.thick_end = 1, target_bed.end

            self.session.add(_Storer(num_id, str(target_bed).encode(), transcript.format("gff3").encode()))
            self.session.commit()


def transfer_by_alignment(ref_pep, target_cdna, target_bed):
    frames = dict()
    # Get the three-frame translation
    for frame in range(3):
        frames[frame] = str(Seq.Seq(str(target_cdna[frame:])).translate(to_stop=False))

    # This will get the best match in the 3-frame translation
    aln, best_frame, score, best_cigar = None, None, float("-inf"), None
    for frame in frames:
        res, cigar = transfer.get_and_prepare_cigar(ref_pep,
                                                    frames[frame],
                                                    open=3,
                                                    extend=1,
                                                    matrix=parasail.blosum85)
        if res.score > score:
            aln, best_frame, score, best_cigar = res, frame, res.score, cigar

    # Now it is time to try to transfer it ... Ignore any deletions at the beginning
    cig_start = 0
    translation_start = 0

    for cig_pos, cig in enumerate(best_cigar):
        le, op = cig
        if not transfer.op_consumes[op][0]:
            # Pass by deletions
            cig_start += 1
            translation_start += best_cigar[cig_start][0]
            continue
        else:
            if transfer.op_consumes[op][1]:
                # anslation_start += best_cigar[cig_start][0]
                break
            else:
                cig_start += 1
                continue

    # This is 0-based; we have to add 1 because we start 1 base after the gap at the beginning
    if translation_start > 0:
        translation_start = 3 * translation_start + 1
    else:
        # We have to account for the frame!
        translation_start = best_frame

    translated = str(Seq.Seq(str(target_cdna[translation_start:])).translate(to_stop=True))

    # Logic to handle when the CDS is broken
    # This is 1-based, so we have to add 1 to
    target_bed.thick_start = translation_start + 1
    end = target_bed.thick_start + len(translated) * 3 - 1

    if translated and translated[0] != ref_pep[0]:
        if translation_start in (0, 1, 2):
            target_bed.phase = translation_start
            target_bed.thick_start = 1
        else:
            target_bed.coding = False
            return target_bed, ()
    elif not translated:
        target_bed.coding = False
        return target_bed, ()

    # Get the coordinates on the original protein

    pep_res, pep_cigar = transfer.get_and_prepare_cigar(ref_pep, translated,
                                                        open=3, extend=1, matrix=parasail.blosum85)

    pep_ref_array, pep_target_array = transfer.create_translation_array(pep_cigar)
    pep_start, pep_end = None, None

    for pos in range(1, len(pep_ref_array) + 1):
        if pep_ref_array[pos - 1] and pep_target_array[pos - 1]:
            if not pep_start:
                pep_start = pos
        pep_end = pos

    # Now check whether we can add the stop codon
    if end + 3 < len(target_cdna):
        end += 3
    else:  # Here we have to presume that it is open.
        end = len(target_cdna)

    # print(translation_start * 3, translated)
    target_bed.thick_end = end
    target_bed.coding = True
    target_bed.transcriptomic = True
    return target_bed, (pep_start, pep_end)


def transfer_cds(transcript, ref_cdna, ref_bed, target_cdna, target_bed):

    if transcript is None:
        return

    transcript.finalize()
    assert target_bed.transcriptomic is True

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

    # We presume it is correct
    pep_coords = (1, len(ref_pep))
    if not replace:
        transcript.attributes["original_cds"] = True
        transcript.attributes["aligner_cds"] = True
        print("CDS transferred correctly for {}".format(ref_bed.chrom))
    else:
        transcript.strip_cds()
        transcript.attributes["original_cds"] = True
        transcript.attributes["aligner_cds"] = False
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
            target_bed, pep_coords = transfer_by_alignment(ref_pep, target_cdna, target_bed)
            # Let's see what happens when we just use the transferred cDNA start

            if target_bed.coding is True and target_bed.invalid is False:
                transcript.attributes["original_pep_coords"] = "{}-{}".format(*pep_coords)
                transcript.attributes["original_pep_complete"] = (pep_coords[0] == 1 and pep_coords[1] == len(ref_pep))
                if (orig_start, orig_end) == (target_bed.thick_start, target_bed.thick_end):
                    print("GMAP corrected the CDS and phase for {}: from {}-{} to {}-{}".format(
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

    if pep_coords:
        pep_coords = (pep_coords[0], pep_coords[1], (pep_coords[0] == 1 and pep_coords[1] == len(ref_pep)))

    return transcript, target_bed, pep_coords


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--bed12", nargs=2, required=True, help="Transcriptomic cDNAs BED12s")
    parser.add_argument("--cdnas", nargs=2, required=True)
    parser.add_argument("-gf", help="GFF3 of the transferred annotation.", required=True)
    parser.add_argument("--out", default=sys.stdout, type=argparse.FileType("wt"))
    parser.add_argument("-ob", "--out-bed", dest="out_bed", required=False,
                        default=None,
                        type=argparse.FileType("wt"))
    parser.add_argument("-p", "--processes", type=int, default=mp.cpu_count())
    args = parser.parse_args()

    cdnas = dict()
    beds = dict()

    gmap_pat = re.compile("\.mrna[0-9]*$")

    for key, cdna, bed in zip(("ref", "target"), args.cdnas, args.bed12):
        cdnas[key] = pyfaidx.Fasta(cdna)
        beds[key] = dict()
        for entry in Bed12Parser(bed, transcriptomic=True):
            if entry.header:
                continue
            beds[key][re.sub(gmap_pat, "", entry.chrom)] = entry

    # Now let us start parsing the GFF3, which we presume being a GMAP GFF3
    transcript = None

    procs = []
    queue = mp.Queue(-1)
    for proc in range(args.processes):
        sq = tempfile.NamedTemporaryFile(mode="wb")
        sq.close()
        sq = sq.name
        _proc = Transferer(sq, queue)
        _proc.start()
        procs.append(_proc)

    # pool = mp.Pool(processes=args.processes)

    tnum = -1
    for line in GFF3(args.gf):
        if line.header is True or line.gene is True:
            print(line, file=args.out)
            continue
        elif line.is_transcript is True:
            if transcript:
                tid = re.sub(gmap_pat, "", transcript.id)
                ref_cdna = str(cdnas["ref"][tid])
                ref_bed = beds["ref"][tid]
                target_cdna = str(cdnas["target"][transcript.id])
                target_bed =  beds["target"][tid]
                tnum += 1
                queue.put((tnum,
                           (transcript,
                            ref_cdna, ref_bed,
                            target_cdna, target_bed)
                           ))
            transcript = Transcript(line)
        elif line.is_exon is True:
            transcript.add_exon(line)

    if transcript:
        tnum += 1
        tid = re.sub(gmap_pat, "", transcript.id)
        ref_cdna = str(cdnas["ref"][tid])
        ref_bed = beds["ref"][tid]
        target_cdna = str(cdnas["target"][transcript.id])
        target_bed = beds["target"][tid]
        queue.put((tnum,
                   (transcript,
                    ref_cdna, ref_bed,
                    target_cdna, target_bed)
                   ))

    queue.put("EXIT")
    [_proc.join() for _proc in procs]

    # Now the printing ...
    # results = dict()

    for proc in procs:
        sq = sqlalchemy.create_engine("sqlite:///{}".format(proc.out_sq))
        for res in sq.execute("select * from storer"):
            num, bed12, gff3 = res
            if args.out_bed is not None:
                print(bed12.decode(), file=args.out_bed)
            print(*gff3.decode().split("\n"), file=args.out, sep="\n")
        os.remove(proc.out_sq)

    return


if __name__ == "__main__":
    main()
