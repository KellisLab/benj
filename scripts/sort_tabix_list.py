#!/usr/bin/env python3

def run(input_files, output, chrom_sizes, buffer_size=100000, threads=8):
    import pysam
    import heapq
    import sys
    FL = [pysam.TabixFile(fname, threads=threads) for fname in input_files]
    contig_list = FL[0].contigs
    for i in range(1, len(FL)):
        if set(contig_list) != set(FL[i].contigs):
            print("Warning: Contigs are not equal", file=sys.stderr)
    for contig in contig_list:
        for start in range(0, chrom_sizes[contig], buffer_size):
            print(contig, start, file=sys.stderr)
            end = start + buffer_size
            heap = []
            for bed in FL:
                for line in bed.fetch(contig, start, end, parser=pysam.asBed()):
                    heapq.heappush(heap, (line.start, str(line)))
            while heap:
                print(heapq.heappop(heap)[1], file=output)

if __name__ == "__main__":
    import argparse
    import sys
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, nargs="+", dest="input_files")
    ap.add_argument("-o", "--output", type=argparse.FileType("w"), default=sys.stdout)
    ap.add_argument("-c", "--chrom-sizes", required=True)
    ap.add_argument("-b", "--buffer-size", type=int, default=100000)
    args = vars(ap.parse_args())
    chrom_sizes = pd.read_csv(args["chrom_sizes"], header=None, sep="\t")
    run(input_files=args["input_files"], output=args["output"], buffer_size=args["buffer_size"], chrom_sizes={k: v for k, v in zip(chrom_sizes[0], chrom_sizes[1])})
