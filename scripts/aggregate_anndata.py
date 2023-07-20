#!/usr/bin/env python3

def add_cell_cycle_info(adata, cell_cycle):
    return 0

def run(metadata, output, directory=[], sample_key="Sample", cell_cycle=None, gtf=None, min_cells_per_sample:int=30, compression:int=9, min_n_genes:int=0, min_total_counts:int=0, use_scrublet:bool=True, qc_vars=[], **kwargs):
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
    scrub_table = {}
    with sw("Reading H5AD files"):
        for sample in tqdm(md.index.values):
            for dname in directory:
                fname = os.path.join(dname, "%s.h5ad" % sample)
                if os.path.isfile(fname):
                    adata = sc.read(fname)
                    if "scrublet" in adata.uns:
                        scrub_table[sample] = adata.uns["scrublet"]
                    for cn in md.columns:
                        adata.obs[cn] = md.loc[sample, cn]
                    if "n_genes_by_counts" in adata.obs.columns and min_n_genes > 0:
                        adata = adata[adata.obs["n_genes_by_counts"] >= min_n_genes, :].copy()
                    if "total_counts" in adata.obs.columns and min_total_counts > 0:
                        adata = adata[adata.obs["total_counts"] >= min_total_counts, :].copy()
                    if adata.shape[0] >= min_cells_per_sample:
                        tbl[sample] = adata
                    else:
                        bad.append(sample)
                    break
            else:
                raise RuntimeError("Sample %s.h5ad did not exist in the directories: %s" % (sample, ",".join(directory)))
    if bad:
        print("Bad samples: ", ",".join(bad))
    total_cells = np.sum([adata.shape[0] for _, adata in tbl.items()])
    with sw("Concatenating %d cells into one AnnData object" % total_cells):
        adata = anndata.concat(tbl, merge="same", uns_merge="same")
        tk = tbl.keys()
        del tbl
        if len(tk) == len(scrub_table.keys()):
            adata.uns["scrublet"] = {"batches": scrub_table,
                                     "batched_by": "Sample"}
    with sw("Calculating statistics"):
        if qc_vars is None:
            qc_vars=[]
        sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True, percent_top=[])
    if "gene_ids" in adata.var.columns:
        with sw("Finding HVGs"):
            try:
                sc.pp.highly_variable_genes(adata, batch_key="Sample", flavor="seurat_v3", n_top_genes=adata.shape[1])
            except ValueError:
                pass
    if use_scrublet and "scrublet" not in adata.uns.keys():
        with sw("Scrublet"):
            sc.external.pp.scrublet(adata, batch_key="Sample")
    with sw("Writing H5AD"):
        adata.write_h5ad(output, compression="gzip", compression_opts=compression)



def read_backed_noX(fname, min_n_genes:int=0, **kwargs):
    """Read an AnnData object except for .X, .raw, and .layers"""
    import h5py
    import anndata
    from anndata.experimental import read_elem
    attrs = ["obs", "obsm", "obsp", "var", "varm", "varp", "uns"]
    with h5py.File(fname, "r") as f:
        adata = anndata.AnnData(**{k: read_elem(f[k]) for k in attrs})
    for cn, val in kwargs.items():
        adata.obs[cn] = val
    if "n_genes_by_counts" in adata.obs.columns and min_n_genes > 0:
        adata = adata[adata.obs["n_genes_by_counts"] >= min_n_genes, :].copy()
    return adata

def run_backed(metadata, output, directory=[], sample_key="Sample", cell_cycle=None, gtf=None, min_cells_per_sample:int=30, compression:int=9, min_n_genes:int=0, use_scrublet:bool=True, qc_vars=[], **kwargs):
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
    with sw("Finding H5AD files"):
        for sample in tqdm(md.index.values):
            for dname in directory:
                fname = os.path.join(dname, "%s.h5ad" % sample)
                if os.path.isfile(fname):
                    import h5py
                    with h5py.File(fname, "r") as f:
                        ashape = f["X"].attrs["shape"]
                    if ashape[0] >= min_cells_per_sample:
                        tbl[sample] = fname
                    else:
                        bad.append(sample)
                    break
            else:
                raise RuntimeError("Sample %s.h5ad did not exist in the directories: %s" % (sample, ",".join(directory)))
    if bad:
        print("Bad samples: ", ",".join(bad))
    with sw("Concatenating H5AD non-X items"):
        adata = anndata.concat({sample: read_backed_noX(fname, min_n_genes=min_n_genes,
                                                        **{k: md.loc[sample, k] for k in md.columns})
                                for sample, fname in tbl.items()},
                               merge="same", uns_merge="same")
    with sw("Writing X-less H5AD to disk"):
        adata.write_h5ad(output, compression="gzip", compression_opts=compression)
    with sw("Loading X (%d, %d)" % adata.shape):
        from anndata.experimental.multi_files import AnnCollection
        adata.X = AnnCollection({k: sc.read(v, backed="r") for k, v in tbl.items()},
                                join_vars="inner")[adata.obs_names, adata.var_names].X
        del tbl
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
    ap.add_argument("--min-cells-per-sample", type=int, default=30, help="Minimum number of cells per sample in order to include")
    ap.add_argument("--use-scrublet", dest="use_scrublet", action="store_true")
    ap.add_argument("--no-use-scrublet", dest="use_scrublet", action="store_false")
    ap.add_argument("--min-n-genes", type=int, default=0)
    ap.add_argument("--min-total-counts", type=int, default=0)
    ap.add_argument("--compression", type=int, default=9)
    ap.add_argument("--backed", dest="backed", action="store_true", help="Use AnnCollection experimental API to load backed anndata objects")
    ap.add_argument("--no-backed", dest="backed", action="store_false")
    ap.set_defaults(scrublet=True, backed=False)
    args = vars(ap.parse_args())
    if args["backed"]:
        run_backed(**args)
    else:
        run(**args)
