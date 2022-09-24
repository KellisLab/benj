#!/usr/bin/env python3
import argparse
import os
import numpy as np
import benj

def gen_umap(adata, batch="Sample", sample_size=0.10, min_components=5, plot=[], resolution=2.):
    from harmonypy import run_harmony
    from umap import UMAP
    import scanpy as sc
    benj.spectral(adata, sample_size=sample_size)
    ns = benj.n_spectral(adata.uns["spectral_eigenvalue"])
    ns = max(ns, min_components)
    print("Using %d components" % ns)
    rep = "X_spectral"
    if batch in adata.obs.columns:
        adj_rep = "%s_harmony" % rep
        adata.obsm[adj_rep] = run_harmony(
            adata.obsm[rep][:, :ns],
            adata.obs,
            batch,
            max_iter_harmony=50
        ).Z_corr.T
        rep = adj_rep
    elif batch is not None:
        print("Batch \"%s\" does not exist.\nUsing uncorrected X_spectral instead." % batch)
    adata.obsm["X_umap"] = UMAP(random_state=0).fit_transform(
        adata.obsm[rep][:, :ns]
    )
    if "Sample" in adata.obs.columns:
        sc.pl.umap(adata, color="Sample", save="_sample.png")
    for item in plot:
        sc.pl.umap(adata, color=item, save="_%s.png" % item)
    sc.pp.neighbors(adata, use_rep=rep, n_pcs=ns)
    sc.tl.leiden(adata, resolution=resolution)
    sc.pl.umap(adata, color="leiden", save="_leiden.png")
    return adata

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-b", "--batch", type=str, default=None)
    ap.add_argument("--min-components", type=int, default=5)
    ap.add_argument("-p", "--plot", nargs="+")
    ap.add_argument("-r", "--resolution", default=2., type=float)
    ap.add_argument("--sample-size", type=float, default=0.10)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    adata = benj.parse_anndata(h5ad=args["input"], **args)
    gen_umap(adata, sample_size=args["sample_size"], batch=args["batch"], min_components=args["min_components"], plot=args["plot"], resolution=args["resolution"])
    adata.write_h5ad(args["output"], compression="gzip")
