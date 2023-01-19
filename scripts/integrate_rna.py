#!/usr/bin/env python3
def integrate(adata, batch=None, hvg=0, use_combat=True, use_harmony=True, use_bbknn=True, plot=None, leiden="overall_clust", resolution=1., dotplot=None, celltypist_model=None, tsv=None, **kwargs):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    if batch not in adata.obs.columns:
        batch = None
    if batch is not None:
        cdf = adata.obs.groupby(batch).count().iloc[:, 0]
        cdf = cdf[cdf >= 3].index.values
        if np.setdiff1d(pd.unique(adata.obs[batch]), cdf) > 0:
            adata = adata[adata.obs[batch].isin(cdf.index), :].copy()
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    adata.raw = adata
    if hvg > 0:
        sc.pp.highly_variable_genes(adata, n_top_genes=hvg, batch_key=batch, subset=True)
    if batch is not None and use_combat:
        sc.pp.combat(adata, batch)
    else:
        sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata)
    if batch is not None and use_harmony:
        sc.external.pp.harmony_integrate(adata, batch, max_iter_harmony=50)
        rep = "X_pca_harmony"
    else:
        rep = "X_pca"
    if batch is not None and use_bbknn:
        sc.external.pp.bbknn(adata, batch_key=batch, use_rep=rep)
    else:
        sc.pp.neighbors(adata, use_rep=rep)
    sc.tl.umap(adata)
    if plot is not None:
        plot = np.union1d(["biosample", "pathology", "log1p_total_counts"], plot)
    else:
        plot = ["biosample", "pathology", "log1p_total_counts"]
    for col in plot:
        if col in adata.obs.columns or col in adata.var.index:
            sc.pl.umap(adata, color=col, save="_%s.png" % col)
    sc.tl.leiden(adata, resolution=resolution, key_added=leiden)
    if tsv is not None:
        adata.obs.loc[:, [leiden]].to_csv(tsv, sep="\t")
    sc.pl.umap(adata, color=leiden, save="_%s_beside.png" % leiden)
    sc.pl.umap(adata, color=leiden, save="_%s_ondata.png" % leiden, legend_loc="on data")
    if dotplot is not None:
        sc.pl.dotplot(adata, var_names=dotplot, groupby=leiden, save="%s.png" % leiden, standard_scale="var")
    for vv in np.intersect1d(["pct_counts_mt", "doublet_score", "log1p_total_counts"], adata.obs.columns):
        sc.pl.violin(adata, vv, groupby=leiden, save="_%s_%s.png" % (leiden, vv))
    sc.tl.dendrogram(adata, groupby=leiden)
    sc.tl.rank_genes_groups(adata, groupby=leiden, method="wilcoxon", pts=True)
    sc.pl.rank_genes_groups_dotplot(adata, save="rgg_%s.png" % leiden)
    sc.pl.rank_genes_groups_matrixplot(adata, save="rgg_%s.png" % leiden)
    sc.pl.rank_genes_groups_heatmap(adata, save="_rgg_%s.png" % leiden)
    adata.write_h5ad(output, compression="gzip")
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
    ap.add_argument("--hvg", default=0, type=int)
    ap.add_argument("--no-use-combat", dest="use_combat", action="store_false")
    ap.add_argument("--use-combat", dest="use_combat", action="store_true")
    ap.add_argument("--no-use-harmony", dest="use_harmony", action="store_false")
    ap.add_argument("--use-harmony", dest="use_harmony", action="store_true")
    ap.add_argument("--no-use-bbknn", dest="use_bbknn", action="store_false")
    ap.add_argument("--use-bbknn", dest="use_bbknn", action="store_true")
    ap.add_argument("--dotplot", nargs="+")
    ap.set_defaults(use_combat=True, use_harmony=True, use_bbknn=True)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    adata = benj.parse_anndata(**args)
    integrate(adata, **args)
