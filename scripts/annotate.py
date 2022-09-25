#!/usr/bin/env python3

import argparse
import os
import benj

def annotate(adata, **kwargs):
    ct = benj.annotate(adata, **kwargs)
    folder = kwargs["folder"]
    if not os.path.isdir(folder):
        os.mkdir(folder)
    ct.to_table(folder, kwargs["prefix"])
    df = ct.predicted_labels
    cluster = kwargs["cluster"]
    df[cluster] = adata.obs[cluster]
    if kwargs["majority_voting"]:
        df = benj.annotate_clusters_from_vote(df, cluster, "majority_voting", kwargs["newcol"])
    else:
        df = benj.annotate_clusters_from_vote(df, cluster, "predicted_labels", kwargs["newcol"])
    df.to_csv(os.path.join(folder, "%s%s.csv.gz" % (kwargs["prefix"], kwargs["newcol"])))
    return df

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-m", "--model")
    ap.add_argument("-f", "--folder", default=os.getcwd())
    ap.add_argument("-p", "--prefix", default="")
    ap.add_argument("--majority-voting", dest="majority_voting", action="store_true")
    ap.add_argument("--no-majority-voting", dest="majority_voting", action="store_false")
    ap.add_argument("--over-clustering", default="")
    ap.add_argument("--cluster", default="leiden")
    ap.add_argument("--newcol", default="CellType")
    ap.set_defaults(majority_voting=True)
    args = benj.parse_args(ap, ["log", "anndata"])
    adata = benj.parse_anndata(h5ad=args["input"], **args)
    annotate(adata=adata, **args)
