#!/usr/bin/env python3

import argparse
import os
import sys

def annotate(adata, model, folder, prefix, majority_voting=True):
    import celltypist
    import numpy as np
    import scanpy as sc
    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)
    ct = celltypist.annotate(adata, model=model, majority_voting=majority_voting)
    ct.to_table(folder, prefix)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-m", "--model", required=True)
    ap.add_argument("-f", "--folder", default=os.getcwd())
    ap.add_argument("-p", "--prefix", default="")
    ap.add_argument("--majority-voting", dest="majority_voting", action="store_true")
    ap.add_argument("--no-majority-voting", dest="majority_voting", action="store_false")
    ap.set_defaults(majority_voting=True)
    args = vars(ap.parse_args())
    if os.path.exists(args["input"]):
        import anndata
        adata = anndata.read(args["input"])
    else:
        print("File \"%s\" does not exist." % args["input"])
        sys.exit(1)
    del args["input"]
    annotate(adata=adata, **args)
