#!/usr/bin/env python3


if __name__ == "__main__":
        import argparse
        import sys, os
        import pandas as pd
        ap = argparse.ArgumentParser()
        ap.add_argument("-i", "--input", dest="csv", nargs="+")
        ap.add_argument("-o", "--output", required=True)
        ap.add_argument("-s", "--sep", default=",")
        args = vars(ap.parse_args())
        out = {}
        for csv in args["csv"]:
                outs_path = os.path.dirname(csv)
                sample_path = os.path.dirname(outs_path)
                out[os.path.basename(sample_path)] = pd.read_csv(csv, thousands=",")
        pd.concat(out).reset_index(level=1, drop=True).to_csv(args["output"], sep=args["sep"])
