#!/usr/bin/env python3

def add_cell_cycle_info(adata, cell_cycle):
    return 0

def run(metadata, output, directory=[], sample_key="Sample", cell_cycle=None, gtf=None, min_cells_per_sample:int=30, compression:int=9, min_n_genes:int=0, use_scrublet:bool=True, qc_vars=[], **kwargs):
    import os
    from tqdm.auto import tqdm
    import pandas as pd
    import numpy as np
    import anndata
    import scanpy as sc
    import benj
    md = pd.read_csv(metadata, sep="\t", index_col=0)
    tbl = {}
    sw = benj.stopwatch()
    bad = []
    if use_scrublet:
        min_n_genes = max(min_n_genes, 3)
    with sw("Reading H5AD files"):
        for sample in tqdm(md.index.values):
            for dname in directory:
                fname = os.path.join(dname, "%s.h5ad" % sample)
                if os.path.isfile(fname):
                    adata = sc.read(fname)
                    for cn in md.columns:
                        adata.obs[cn] = md.loc[sample, cn]
                    if adata.shape[0] >= min_cells_per_sample:
                        if "n_genes_by_counts" in adata.obs.columns and min_n_genes > 0:
                            adata = adata[adata.obs["n_genes_by_counts"] >= min_n_genes, :].copy()
                        tbl[sample] = adata
                    else:
                        bad.append(sample)
                    break
            else:
                raise RuntimeError("Sample %s.h5ad did not exist in the directories: %s" % (sample, ",".join(directory)))
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
    ap.add_argument("-d", "--directory", nargs="+", help="Directories with {Sample}.h5ad files inside", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--qc-vars", nargs="+")
    ap.add_argument("--min-cells-per-sample", type=int, default=30)
    ap.add_argument("--use-scrublet", dest="use_scrublet", action="store_true")
    ap.add_argument("--no-use-scrublet", dest="use_scrublet", action="store_false")
    ap.add_argument("--min-n-genes", type=int, default=0)
    ap.add_argument("--compression", type=int, default=9)
    ap.set_defaults(scrublet=True)
    args = vars(ap.parse_args())
    run(**args)
