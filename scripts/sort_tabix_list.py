#!/usr/bin/env python3
import functools

@functools.total_ordering
class BedLine:
    def __init__(self, line, **kwargs):
        self.line = line.strip()
        S = self.line.split("\t")
        self.chrom = S[0]
        self.start = int(S[1])
        self.end = int(S[2])
        self.kwargs = kwargs
    def __repr__(self):
        return "BedLine(chrom=%s start=%d end=%d)" % (self.chrom, self.start, self.end)
    def __lt__(self, other):
        return self.chrom < other.chrom or self.start < other.start or self.end < other.end
    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

def run(input_files, output, buffer_size=10):
    import gzip
    import heapq
    FL = [gzip.open(fname, "r") for fname in input_files]
    pq = []
    for index, fname in enumerate(FL):
        for _ in range(buffer_size):
            next_line = fname.readline().decode("utf-8")
            while next_line and next_line.startswith("#"):
                next_line = fname.readline().decode("utf-8")
            if next_line:
                heapq.heappush(pq, BedLine(next_line, index=index))
    while heapq:
        item = heapq.heappop(pq)
        index = item.kwargs["index"]
        print(item.line, file=output)
        if not FL[index].closed:
            next_line = FL[index].readline().decode("utf-8")
            if next_line:
                heapq.heappush(pq, BedLine(next_line, index=index))
            else:
                FL[index].close()

if __name__ == "__main__":
    import argparse
    import sys
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, nargs="+", dest="input_files")
    ap.add_argument("-o", "--output", type=argparse.FileType("w"), default=sys.stdout)
    ap.add_argument("-b", "--buffer-size", type=int, default=10)
    args = vars(ap.parse_args())
    run(input_files=args["input_files"], output=args["output"])
