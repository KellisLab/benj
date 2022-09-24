#!/usr/bin/env python3

import argparse
import os
import benj

def annotate(adata, **kwargs):
    import celltypist
    import numpy as np
    import scanpy as sc
    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)
    folder = kwargs["folder"]
    valid_k = ["over_clustering", "model", "majority_voting"]
    ctargs = {k: v for k, v in kwargs.items() if k in valid_k and v is not None and v != ""}
    ct = celltypist.annotate(adata, **ctargs)
    if not os.path.isdir(folder):
        os.mkdir(folder)
    ct.to_table(folder, kwargs["prefix"])

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-m", "--model")
    ap.add_argument("-f", "--folder", default=os.getcwd())
    ap.add_argument("-p", "--prefix", default="")
    ap.add_argument("--majority-voting", dest="majority_voting", action="store_true")
    ap.add_argument("--no-majority-voting", dest="majority_voting", action="store_false")
    ap.add_argument("--over-clustering", default="")
    ap.set_defaults(majority_voting=True)
    args = benj.parse_args(ap, ["log", "anndata"])
    adata = benj.parse_anndata(h5ad=args["input"], **args)
    annotate(adata=adata, **args)
