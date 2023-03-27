#!/usr/bin/env python3

def integrate(adata, output=None, batch=None, use_harmony=True, use_bbknn=False,
              leiden="overall_clust", resolution=1., prefix="",
              tsv=None, min_dist:float=0.3, compression:int=9,
              min_n_cells_by_counts:int=2, cor_cutoff:float=0.8,
              plot=[],
              qc_cols=["log1p_total_counts"], sw=None, **kwargs):
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import muon as mu
    from muon import atac as ac
    import benj
    if sw is None:
        sw = benj.stopwatch()
    if min_n_cells_by_counts > 0 and "n_cells_by_counts" in adata.var.columns:
        with sw("Filtering cells by counts >= %d" % min_n_cells_by_counts):
            mu.pp.filter_var(adata, "n_cells_by_counts", lambda x: x >= min_n_cells_by_counts)
    if batch not in adata.obs.columns:
        batch = None
    if batch is not None:
        ### Filter batch
        cdf = adata.obs.groupby(batch).count().iloc[:, 0]
        cdf = cdf[cdf >= 3].index.values
        if len(np.setdiff1d(pd.unique(adata.obs[batch]), cdf)) > 0:
            mu.pp.filter_obs(adata, batch, lambda x: x.isin(list(cdf)))
        if len(pd.unique(adata.obs[batch])) <= 1:
            batch = None
        else:
            adata.obs[batch] = adata.obs[batch].values.astype(str)
    if "raw" not in adata.layers:
        with sw("Copying .X to .layers[\"raw\"]"):
            adata.layers["raw"] = adata.X.copy()
    with sw("Re-calculating qc metrics"):
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=True, inplace=True, layer="raw")
    with sw("Running TF-IDF"):
        adata.X = adata.layers["raw"].copy()
        ac.pp.tfidf(adata)
    with sw("Running LSI"):
        ac.tl.lsi(adata)
        adata.X = adata.X.astype(np.float32)
    for col in qc_cols:
        with sw("Correlating column \"%s\" with LSI" % col):
            from scipy.stats import pearsonr
            cor = np.zeros(len(adata.uns["lsi"]["stdev"]))
            for i in range(len(cor)):
                cor[i], _ = pearsonr(adata.obsm["X_lsi"][:, i], adata.obs[col])
            cor_flag = np.abs(cor) < cor_cutoff
            if np.sum(cor_flag) > 0:
                print("Removing components:", np.ravel(np.where(~cor_flag)))
                print("  with correlations", cor[~cor_flag])
            adata.obsm["X_lsi"] = adata.obsm["X_lsi"][:, cor_flag].astype(np.float32)
            adata.varm["LSI"] = adata.varm["LSI"][:, cor_flag].astype(np.float32)
            adata.uns["lsi"]["stdev"] = adata.uns["lsi"]["stdev"][cor_flag]
    use_rep="X_lsi"
    n_pcs=len(adata.uns["lsi"]["stdev"])
    if batch is not None and use_harmony:
        with sw("Running Harmony"):
            use_rep_adj = "%s_harmony" % use_rep
            sc.external.pp.harmony_integrate(adata, batch, basis=use_rep, adjusted_basis=use_rep_adj)
            use_rep = use_rep_adj
    if batch is not None and use_bbknn:
        with sw("Running BB-KNN"):
            sc.external.pp.bbknn(adata, batch, use_rep=use_rep, n_pcs=n_pcs)
    else:
        with sw("Running nearest neighbors"):
            sc.pp.neighbors(adata, use_rep=use_rep, n_pcs=n_pcs)
    with sw("Computing UMAP"):
        sc.tl.umap(adata, min_dist=min_dist)
    with sw("Plotting"):
        if plot is not None:
            plot = np.union1d(plot, ["log1p_total_counts", "TSSEnrichment", "tss_score"])
        else:
            plot = ["log1p_total_counts", "TSSEnrichment", "tss_score"]
        for col in plot:
            if col in adata.obs.columns:
                sc.pl.umap(adata, color=col, save="_%s.png" % col)
            elif col in adata.var.index:
                ac.pl.umap(adata, color=col, save="_%s.png" % col, use_raw=False)
    with sw("Clustering cells"):
        sc.tl.leiden(adata, key_added=leiden, resolution=resolution)
        adata.obs[leiden] = ["%s%s" % (prefix, v) for v in adata.obs[leiden].values.astype(str)]
    if tsv is not None:
        with sw("Writing TSV"):
            cols = np.intersect1d([leiden, "majority_voting", "predicted_labels"], adata.obs.columns)
            adata.obs.loc[:, cols].to_csv(tsv, sep="\t")
    with sw("Plotting %s" % leiden):
        sc.pl.umap(adata, color=leiden, save="_%s_beside.png" % leiden)
        sc.pl.umap(adata, color=leiden, save="_%s_ondata.png" % leiden, legend_loc="on data", legend_fontsize=4)
    with sw("Writing to H5AD"):
        adata.write_h5ad(output, compression="gzip", compression_opts=compression)
    return adata

if __name__ == "__main__":
    import benj
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", dest="h5ad", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-t", "--tsv")
    ap.add_argument("-b", "--batch", type=str, default=None)
    ap.add_argument("-p", "--plot", nargs="+")
    ap.add_argument("-r", "--resolution", default=1., type=float)
    ap.add_argument("-l", "--leiden", default="overall_clust")
    ap.add_argument("--prefix", default="C")
    ap.add_argument("--cor-cutoff", type=float, default=0.8)
    ap.add_argument("--no-use-harmony", dest="use_harmony", action="store_false")
    ap.add_argument("--use-harmony", dest="use_harmony", action="store_true")
    ap.add_argument("--no-use-bbknn", dest="use_bbknn", action="store_false")
    ap.add_argument("--use-bbknn", dest="use_bbknn", action="store_true")
    ap.add_argument("--qc-cols", nargs="+", default=["log1p_total_counts"])
    ap.add_argument("--min-dist", type=float, default=0.5)
    ap.add_argument("--compression", type=int, default=9)
    ap.set_defaults(use_harmony=True, use_bbknn=False)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    adata = benj.parse_anndata(**args)
    integrate(adata, **args)
