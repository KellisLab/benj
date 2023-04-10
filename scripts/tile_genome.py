#!/usr/bin/env python3

def run(chrNameLength, output, tile:int=500, **kwargs):
    for rline in chrNameLength:
        S = rline.strip().split("\t")
        chrom = S[0]
        length = int(S[1])
        for begin in range(0, length, tile):
            end = min(begin + tile, length)
            begin = max(1, begin)
            print("%s\t%d\t%d" % (chrom, begin, end), file=output)

if __name__ == "__main__":
    import argparse
    import sys
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--chrNameLength", type=argparse.FileType("r"), required=True)
    ap.add_argument("-t", "--tile", type=int, default=500)
    ap.add_argument("-o", "--output", type=argparse.FileType("w"), default=sys.stdout)
    args = vars(ap.parse_args())
    run(**args)
