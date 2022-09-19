#!/usr/bin/env python3

import argparse

def cr2h5ad(infile, outfile, gex_only=False):
    import scanpy as sc
    import numpy as np
    import scipy.sparse
    adata = sc.read_10x_h5(infile, gex_only=gex_only)
    dt = adata.X
    if not np.issubtype(dt, np.integer):
        dat = adata.X
        if scipy.sparse.issparse(adata.X):
            dat = adata.X.data
        if np.all(np.mod(dat, 1) == 0):
            adata.X = adata.X.astype(np.uint32)
    adata.var_names_make_unique()
    adata.write_h5ad(outfile, compression="gzip")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--gex-only", dest="gex_only", action="store_true")
    ap.set_defaults(gex_only=False)
    args = vars(ap.parse_args())
    cr2h5ad(infile=args["input"], outfile=args["output"], gex_only=args["gex_only"])
