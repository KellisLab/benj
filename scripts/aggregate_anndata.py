#!/usr/bin/env python3

def add_cell_cycle_info(adata, cell_cycle):
    return 0

def run(metadata, output, directory=".", sample_key="Sample", cell_cycle=None, gtf=None, min_cells_per_sample:int=30, compression:int=9, use_scrublet:bool=True, qc_vars=[], **kwargs):
    import os
    from tqdm.auto import tqdm
    import pandas as pd
    import numpy as np
    import scanpy as sc
    import benj
    md = pd.read_csv(metadata, sep="\t", index_col=0)
    tbl = []
    sw = benj.stopwatch()
    bad = []
    with sw("Reading H5AD files"):
        for sample in tqdm(md.index.values):
            fname = os.path.join(directory, "%s.h5ad" % sample)
            if os.path.isfile(fname):
                adata = sc.read(fname)
                for cn in md.columns:
                    adata.obs[cn] = md.loc[sample, cn]
                if adata.shape[0] >= min_cells_per_sample:
                    tbl[sample] = adata
                else:
                    bad.append(sample)
    if bad:
        print("Bad samples: ", ",".join(bad))
    with sw("Concatenating AnnData objects"):
        adata = anndata.concat(tbl, merge="same", uns_merge="same")
    with sw("Calculating statistics"):
        if qc_vars is None:
            qc_vars=[]
        sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True, percent_top=[])
    if use_scrublet:
        with sw("Scrublet"):
            sc.external.pp.scrublet(adata, batch_key="Sample")
    with sw("Writing H5AD"):
        adata.write_h5ad(output, compression="gzip", compression_opts=compression)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-m", "--metadata", required=True, help="Pre-formatted metadata TSV with index as sample names")
    ap.add_argument("-d", "--directory", default=".", help="Directory with {Sample}.h5ad files inside")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--qc-vars", nargs="+")
    ap.add_argument("--min-cells-per-sample", type=int, default=30)
    ap.add_argument("--use-scrublet", dest="scrublet", action="store_true")
    ap.add_argument("--no-use-scrublet", dest="scrublet", action="store_false")
    ap.add_argument("--compression", type=int, default=9)
    ap.set_defaults(scrublet=True)
    args = vars(ap.parse_args())
    run(**args)
