#!/usr/bin/env python3


import sys
import argparse
from eiliftover.util import transfer_utilities as transfer
import re
import pyfaidx
from Mikado.parsers.bed12 import BED12, Bed12Parser
from Mikado.parsers import to_gff
from Mikado.transcripts import Transcript
from Mikado.loci import Gene
from Mikado.utilities.log_utils import create_default_logger, create_queue_logger, create_null_logger
from Bio import Seq
import parasail
import multiprocessing as mp
import logging
import logging.handlers
import sqlalchemy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Integer, Column, BLOB
from sqlalchemy.orm import sessionmaker
import tempfile
import os
from collections import defaultdict
import operator


__doc__ = """Script to try to translate the CDS from one coordinate system to another."""

transfer_base = declarative_base()

logging_queue = mp.JoinableQueue(-1)
log_queue_handler = logging.handlers.QueueHandler(logging_queue)
log_queue_handler.setLevel(logging.DEBUG)


class _Storer(transfer_base):

    __tablename__ = "storer"

    id = Column(Integer, primary_key=True)
    bed = Column(BLOB)
    gff3 = Column(BLOB)

    def __init__(self, id, bed, gff3):

        self.id, self.bed, self.gff3 = id, bed, gff3


class Transferer(mp.Process):

    def __init__(self, out_sq, queue, verbosity="INFO"):

        super().__init__()
        self.out_sq = out_sq
        self.logging_queue = logging_queue
        self.logger = create_default_logger("")
        self.log_level = verbosity
        create_queue_logger(self)
        self.engine = sqlalchemy.create_engine("sqlite:///{}".format(out_sq))
        transfer_base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        self.queue = queue

    def run(self):

        self.logger.debug("Starting to wait in the queue")
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
            transcript, target_bed, pep_coords = transfer_cds(transcript, ref_cdna, ref_bed, target_cdna, target_bed,
                                                              logger=self.logger)
            target_bed.name = "ID={};coding={}".format(target_bed.name, target_bed.coding)
            if target_bed.coding:
                target_bed.name = target_bed.name + \
                                  ";phase={};original_pep_coords={}-{};original_pep_complete={}".format(
                                      target_bed.phase, *pep_coords)
            else:
                target_bed.thick_start, target_bed.thick_end = 1, target_bed.end
            gene = Gene(transcript)
            gene.finalize()

            self.session.add(_Storer(num_id, str(target_bed).encode(), gene.format("gff3").encode()))
            self.session.commit()


def transfer_by_alignment(ref_pep, target_cdna, target_bed, logger=create_null_logger()):
    frames = dict()
    # Get the three-frame translation
    logger.debug("Phase for %s: %s", target_bed.name, target_bed.phase)
    for frame in range(3):
        frames[frame] = str(Seq.Seq(str(target_cdna[frame:])).translate(to_stop=False))

    # This will get the best match in the 3-frame translation
    frame_res = dict()
    for frame in frames:
        res, cigar = transfer.get_and_prepare_cigar(ref_pep,
                                                    frames[frame],
                                                    open=3,
                                                    extend=1,
                                                    matrix=parasail.blosum85)
        frame_res[frame] = (res, cigar)
    # Now it is time to try to transfer it ... Ignore any deletions at the beginning
    cig_start = 0
    translation_start = 0
    logger.debug("Frames for %s (phase %s): %s", target_bed.name, target_bed.phase, frame_res)
    best_frame = sorted(frame_res.keys(), key=lambda k: frame_res[k][0].score, reverse=True)[0]
    best_cigar = frame_res[best_frame][1]
    logger.debug("Best frame for %s: %s (cigar: %s)", target_bed.name, best_frame, best_cigar)

    for cig_pos, cig in enumerate(best_cigar):
        le, op = cig
        if not transfer.op_consumes[op][0]:
            # Pass by deletions
            translation_start += best_cigar[cig_start][0]
            cig_start += 1
            continue
        else:
            if transfer.op_consumes[op][1]:
                # anslation_start += best_cigar[cig_start][0]
                break
            else:
                cig_start += 1
                continue

    # This is 0-based; we have to add 1 because we start 1 base after the gap at the beginning
    logger.debug("Translation start for %s: %s; phase: %s", target_bed.name, translation_start, target_bed.phase)
    if translation_start > 0:
        translation_start = 3 * translation_start + best_frame
    else:
        # We have to account for the frame!
        translation_start = best_frame

    translated = str(Seq.Seq(str(target_cdna[translation_start:])).translate(to_stop=(ref_pep.count("*") <= 1)))

    # Logic to handle when the CDS is broken
    # This is 1-based, so we have to add 1 to
    target_bed.thick_start = translation_start + 1
    end = target_bed.thick_start + len(translated) * 3 - 1

    logger.debug("Phase for %s: %s", target_bed.name, target_bed.phase)
    if translated and translated[0] != ref_pep[0]:
        if translation_start in (0, 1, 2):
            target_bed.phase = translation_start
            target_bed.thick_start = 1
        else:
            target_bed.coding = False
            return target_bed, (None, None, False)
    elif not translated:
        target_bed.coding = False
        return target_bed, (None, None, False)

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
    logger.debug("Phase for %s: %s", target_bed.name, target_bed.phase)
    return target_bed, (pep_start, pep_end)


def transfer_cds(transcript: Transcript, ref_cdna: str, ref_bed: BED12,
                 target_cdna: str, target_bed: BED12, logger=create_null_logger()):

    if transcript is None:
        return transcript, target_bed, (None, None, False)

    transcript.finalize()
    assert target_bed.transcriptomic is True

    logger.debug("Starting with %s, phases: %s (BED %s)", transcript.id, transcript.phases, target_bed.phase)

    if ref_bed.coding is False:
        logger.debug("%s is non coding, returning immediately.", transcript.id, transcript.phases)
        transcript.attributes["aligner_cds"] = False
        transcript.attributes["was_coding"] = transcript.is_coding
        target_bed.coding = False
        transcript.strip_cds()
        pep_coords = (None, None, True)
    else:
        original_start, original_end = target_bed.thick_start, target_bed.thick_end
        original_phase, original_phases = target_bed.phase, transcript.phases.copy()
        ref_pep = str(Seq.Seq(str(ref_cdna[ref_bed.thick_start - 1:ref_bed.thick_end])).translate(to_stop=False))

        ref_has_multiple_stops = False
        if ref_pep.count("*") == 0:
            pass
        elif abs(ref_pep.index("*") * 3 - ref_bed.cds_len) in (0, 3):
            ref_pep = ref_pep[:ref_pep.index("*")]  # This is the "good" case: the CDS is correct.
        else:
            ref_has_multiple_stops = True
            logger.warning(
                "The sequence of %s has in frame stop codons. Adjusting the program to take this into account.",
                ref_bed.name)

        logger.debug("%s now has phases: %s (%s)", transcript.id, transcript.phases, target_bed.phase)
        target_bed, pep_coords = transfer_by_alignment(ref_pep, target_cdna, target_bed, logger=logger)
        logger.debug("%s now has phases: %s; target bed: %s", transcript.id, transcript.phases, target_bed.phase)
        pep_coords = (pep_coords[0], pep_coords[1], (pep_coords[0] == 1 and pep_coords[1] == len(ref_pep)))

        if target_bed.thick_start == original_start and target_bed.thick_end == original_end:
            transcript.attributes["aligner_cds"] = True
            logger.debug("%s now has phases: %s", transcript.id, transcript.phases)
        else:
            transcript.attributes["aligner_cds"] = False
            transcript.strip_cds()
            if target_bed.coding is True:
                transcript.load_orfs([target_bed])

        logger.debug("%s now has phases: %s", transcript.id, transcript.phases)
        # Now we have to decide whether the transcript has the "original" CDS or not
        result, cigar = transfer.get_and_prepare_cigar(str(ref_cdna), str(target_cdna))
        ref_array, target_array = transfer.create_translation_array(cigar)
        try:
            target_start = target_array[ref_array.index(ref_bed.thick_start)]
        except IndexError:
            target_start = target_bed.start
        try:
            target_end = target_array[ref_array.index(ref_bed.thick_end)]
        except IndexError:
            target_end = target_bed.end

        if target_start == target_bed.thick_start and target_end == target_bed.thick_end:
            transcript.attributes["original_cds"] = True
        else:
            transcript.attributes["original_cds"] = False

        if ref_cdna == target_cdna:
            logger.debug("%s now has phases: %s", transcript.id, transcript.phases)
            if transcript.is_coding is False:
                raise AssertionError("{} not coding".format(transcript.id))
            elif transcript.attributes["original_cds"] is False:
                raise AssertionError("\n".join([str(_) for _ in [transcript.id,
                    (target_bed.thick_start, target_start, target_bed.thick_start == target_start),
                                            (target_bed.thick_end, target_end, target_bed.thick_end == target_end),
                                            target_bed.thick_start == target_start and target_bed.thick_end == target_end
                                            ]]))

    return transcript, target_bed, pep_coords


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--bed12", nargs=2, required=True, help="Transcriptomic cDNAs BED12s")
    parser.add_argument("--cdnas", nargs=2, required=True)
    parser.add_argument("-gf", help="GFF3/BED12 of the transferred annotation.", required=True)
    parser.add_argument("--out", default=sys.stdout, type=argparse.FileType("wt"))
    parser.add_argument("-ob", "--out-bed", dest="out_bed", required=False,
                        default=None,
                        type=argparse.FileType("wt"))
    log = parser.add_mutually_exclusive_group()
    log.add_argument("-q", "--quiet", default=False, action="store_true")
    log.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument("-p", "--processes", type=int, default=mp.cpu_count())
    args = parser.parse_args()

    logger = create_default_logger("master")
    verbosity = "INFO"
    if args.verbose is True:
        verbosity = "DEBUG"
    elif args.quiet is True:
        verbosity = "WARNING"

    listener = logging.handlers.QueueListener(logging_queue, logger)
    listener.propagate = False
    listener.start()
    logger.setLevel(verbosity)

    cdnas = dict()
    beds = dict()
    beds["ref"] = dict()
    beds["target"] = dict()

    gmap_pat = re.compile("\.mrna[0-9]*$")

    logger.info("Loading reference cDNAS")
    cdnas["ref"] = pyfaidx.Fasta(args.cdnas[0])
    logger.info("Loading target cDNAS")
    cdnas["target"] = pyfaidx.Fasta(args.cdnas[1])
    logger.info("Loaded cDNAs")
    logger.info("Loading reference BED12")
    for entry in Bed12Parser(args.bed12[0], transcriptomic=True):
        if entry.header:
            continue
        name = entry.chrom
        if name in beds["ref"]:
            raise KeyError("Duplicated ID for the reference: {}".format(name))
        if name not in cdnas["ref"]:
            raise KeyError("Reference {} not found in the cDNAs!".format(name))
        beds["ref"][name] = entry

    logger.info("Loading target BED12")
    beds["target"] = defaultdict(dict)
    for entry in Bed12Parser(args.bed12[1], transcriptomic=True):
        # Now, here we have to account for the fact that there *might* be multiple alignments
        name = re.sub(gmap_pat, "", entry.chrom)
        if entry.chrom not in cdnas["target"]:
            raise KeyError("Target {} not found in the cDNAs!".format(entry.chrom))
        beds["target"][name][entry.chrom] = entry
    logger.info("Loaded BED12s")

    # Now let us start parsing the GFF3, which we presume being a GMAP GFF3
    transcript = None

    logger.info("Launching sub-processes")
    procs = []
    queue = mp.Queue(-1)
    for proc in range(args.processes):
        sq = tempfile.NamedTemporaryFile(mode="wb")
        sq.close()
        sq = sq.name
        _proc = Transferer(sq, queue, verbosity=verbosity)
        _proc.start()
        procs.append(_proc)
    logger.info("Launched sub-processes, starting parsing annotation")

    # pool = mp.Pool(processes=args.processes)

    tnum = -1
    if args.gf.endswith(("bed12", "bed")):
        parser = Bed12Parser(args.gf, transcriptomic=False)
        for line in parser:
            if line.header:
                continue
            else:
                transcript = Transcript(line)
                tid = re.sub(gmap_pat, "", transcript.id)
                logger.debug("Found %s", tid)
                ref_cdna = str(cdnas["ref"][tid])
                ref_bed = beds["ref"][tid]
                target_cdna = str(cdnas["target"][transcript.id])
                target_bed = beds["target"][tid][transcript.id]
                tnum += 1
                logger.debug("Submitting %s", tid)
                queue.put((tnum,
                           (transcript,
                            ref_cdna, ref_bed,
                            target_cdna, target_bed)
                           ))
            if tnum >= 10 ** 4 and tnum % 10 ** 4 == 0:
                logger.info("Parsed {} transcripts", tnum)
        logger.info("Finished parsing input genomic BED file")
    else:
        parser = to_gff(args.gf)

        for pos, line in enumerate(parser):
            if line.header is True:  # or (not isinstance(line, BED12) and line.is_gene is True):
                if str(line) == "###":
                    continue
                try:
                    print(line, file=args.out)
                except IndexError:
                    raise IndexError(line._line)
                continue
            elif not isinstance(line, BED12) and line.is_gene is True:
                continue
            elif line.is_transcript is True:
                if transcript:
                    if transcript.alias is None:
                        tid = re.sub(gmap_pat, "", transcript.id)
                    else:
                        tid = re.sub(gmap_pat, "", transcript.alias)
                    ref_cdna = str(cdnas["ref"][tid])
                    ref_bed = beds["ref"][tid]
                    target_cdna = str(cdnas["target"][transcript.id])
                    store = beds["target"].get(tid, None)
                    if store is None:
                        raise KeyError((tid, beds["target"].keys()))
                    target_bed = store.get(transcript.id, None)
                    if target_bed is None:
                        raise KeyError((tid, store.keys()))
                    tnum += 1
                    queue.put((tnum,
                               (transcript,
                                ref_cdna, ref_bed,
                                target_cdna, target_bed)
                               ))
                try:
                    transcript = Transcript(line)
                except (ValueError, TypeError):
                    raise ValueError((pos, line))
            elif line.is_exon is True:
                transcript.add_exon(line)
            if tnum >= 10 ** 4 and tnum % 10 ** 4 == 0:
                logger.info("Parsed {} transcripts", tnum)

        if transcript:
            tnum += 1
            tid = re.sub(gmap_pat, "", transcript.id)
            ref_cdna = str(cdnas["ref"][tid])
            ref_bed = beds["ref"][tid]
            target_cdna = str(cdnas["target"][transcript.id])
            target_bed = beds["target"][tid][transcript.id]
            queue.put((tnum,
                       (transcript,
                        ref_cdna, ref_bed,
                        target_cdna, target_bed)
                       ))
        logger.info("Finished parsing input genomic GF file")

    queue.put("EXIT")
    logger.info("Waiting for subprocesses to finish")
    [_proc.join() for _proc in procs]

    # Now the printing ...
    # results = dict()

    logger.info("Subprocesses finished, printing")
    for proc in procs:
        sq = sqlalchemy.create_engine("sqlite:///{}".format(proc.out_sq))
        for res in sq.execute("select * from storer"):
            num, bed12, gff3 = res
            if args.out_bed is not None:
                print(bed12.decode(), file=args.out_bed)
            print(*gff3.decode().split("\n"), file=args.out, sep="\n")
        os.remove(proc.out_sq)

    logger.info("Finished!")
    return


if __name__ == "__main__":
    main()
