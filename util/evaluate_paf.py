from eiliftover.alignments.cminimap import Stitcher, MiniAlignment
import argparse
import sys
import gzip


def main():
    
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-f", "--flank", default=0, type=int)
    parser.add_argument("minimap")
    parser.add_argument("out", nargs="?", default=sys.stdout, type=argparse.FileType("wt"))
    args = parser.parse_args()

    stitch = Stitcher(flank=args.flank)

    if args.minimap.endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    with opener(args.minimap, mode="rt") as mini:
        for line in mini:
            line = MiniAlignment(line)
            if line.query == stitch.query or stitch.query is None:
                stitch.add(line)
            else:
                # print("Finding the best complete alignment for", stitch.query)
                complete = stitch.best_complete_alignment
                if complete is not None:
                    # print("Found a complete alignment for", stitch.query)
                    print(complete.location.transcript)
                    # print(complete.location["transcript"].decode("utf-8"))

                stitch = Stitcher(flank=args.flank)
                stitch.add(line)

    complete = stitch.best_complete_alignment
    if complete is not None:
        print(complete.location.transcript)

main()


