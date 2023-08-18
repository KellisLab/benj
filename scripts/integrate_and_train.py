#!/usr/bin/env python3

if __name__ == "__main__":
        import argparse
        import os
        import benj
        ap = argparse.ArgumentParser()
        ap.add_argument("-i", "--input", dest="h5ad", required=True)
        ap.add_argument("-o", "--output", required=True)
        ap.add_argument("--labels", required=True, nargs="+")
        ap.add_argument("-b", "--batch", required=True)
        ap.add_argument("--hvg", type=int, default=5000)
        ap.add_argument("--compression", type=int, default=6)
        ap.add_argument("--with-mean", dest="with_mean", action="store_true")
        ap.add_argument("--without-mean", dest="with_mean", action="store_false")
        ap.set_defaults(with_mean=False)
        args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
        if "subset" not in args or args["subset"] is None:
                args["subset"] = []
        sw = benj.stopwatch()
        with sw("Reading H5AD"):
                adata = benj.parse_anndata(**args)
        adata = benj.integrate_rna(adata,
                                   batch=args["batch"], hvg=args["hvg"], use_scaling=True, use_harmony=True, use_bbknn=False, use_rgg=False, plot=args["labels"], target_sum=1e4,
                                   output=args["output"], compression=args["compression"])
        with sw("Training celltypist for " + ",".join(label)):
                import celltypist
                import scanpy as sc
                ct = celltypist.train(adata.raw.to_adata(), labels=adata.obs.loc[:, labels], genes=adata.var_names, n_jobs=-1, with_mean=args["with_mean"])
                ct.write(os.path.join(sc.settings.figdir, f"celltypist_{label}.pkl"))
