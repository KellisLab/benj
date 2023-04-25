#!/usr/bin/env python3
def integrate(adata, output=None, batch=None, hvg=0, use_combat=False, use_scaling=False, use_harmony=True, use_bbknn=True, plot=None, leiden="overall_clust", resolution=1., min_dist:float=0.5, dotplot=None, celltypist=None, tsv=None, rgg_ng=5, max_iter_harmony:int=50, prefix="", sw=None, use_rgg:bool=True, target_sum:int=None, compression:int=9, **kwargs):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import benj
    if sw is None:
        sw = benj.stopwatch()
    if batch not in adata.obs.columns:
        batch = None
    if batch is not None:
        cdf = adata.obs.groupby(batch).count().iloc[:, 0]
        cdf = cdf[cdf >= 3].index.values
        if len(np.setdiff1d(pd.unique(adata.obs[batch]), cdf)) > 0:
            adata = adata[adata.obs[batch].isin(cdf), :].copy()
        if len(pd.unique(adata.obs[batch])) == 0:
            batch = None
    print("Working with %d cells" % adata.shape[0])
    if "raw" in adata.layers:
        with sw("Copying .layers[\"raw\"] to .X"):
            adata.X = adata.layers["raw"].copy()
    else:
        with sw("Copying .X to .layers[\"raw\"]"):
            adata.layers["raw"] = adata.X.copy()
    with sw("Normalizing data"):
        if target_sum is not None and target_sum > 0:
            sc.pp.normalize_total(adata, target_sum=target_sum)
        else:
            print("Using median normalization")
            sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.raw = adata
    if hvg > 0:
        with sw("Calculating %d HVG" % hvg):
            sc.pp.highly_variable_genes(adata, n_top_genes=hvg, batch_key=batch, subset=True)
    if batch is not None and use_combat:
        with sw("Running ComBat"):
            sc.pp.combat(adata, batch)
    elif use_scaling:
        with sw("Scaling data"):
            sc.pp.scale(adata, max_value=10)
    with sw("Running PCA"):
        sc.pp.pca(adata)
    if batch is not None and use_harmony:
        with sw("Running Harmony"):
            sc.external.pp.harmony_integrate(adata, batch, max_iter_harmony=max_iter_harmony)
            rep = "X_pca_harmony"
    else:
        rep = "X_pca"
    if batch is not None and use_bbknn:
        with sw("Running BBKNN"):
            sc.external.pp.bbknn(adata, batch_key=batch, use_rep=rep)
    else:
        with sw("Running neighbors"):
            sc.pp.neighbors(adata, use_rep=rep)
    with sw("Running UMAP"):
        sc.tl.umap(adata, min_dist=min_dist)
    if plot is not None:
        plot = np.union1d(["biosample", "pathology", "log1p_total_counts"], plot)
    else:
        plot = ["biosample", "pathology", "log1p_total_counts"]
    for col in plot:
        if col in adata.obs.columns or col in adata.var.index:
            sc.pl.umap(adata, color=col, save="_%s.png" % col)
    with sw("Running Leiden"):
        sc.tl.leiden(adata, resolution=resolution, key_added=leiden)
        adata.obs[leiden] = ["%s%s" % (prefix, v) for v in adata.obs[leiden].values.astype(str)]
    if celltypist is not None:
        with sw("Annotating from CellTypist"):
            from celltypist import annotate
            if target_sum == 10000:
                ct = annotate(adata.raw.to_adata(), majority_voting=True, over_clustering=leiden, model=celltypist)
            else:
                xdata = anndata.AnnData(adata.layers["raw"], obs=adata.obs, var=adata.var, obsp=adata.obsp)
                sc.pp.normalize_total(xdata, target_sum=10000)
                sc.pp.log1p(xdata)
                ct = annotate(xdata, majority_voting=True, over_clustering=leiden, model=celltypist)
                del xdata
            for cn in ct.predicted_labels.columns:
                adata.obs[cn] = ct.predicted_labels[cn]
            del ct
    if tsv is not None:
        cols = np.intersect1d([leiden, "majority_voting", "predicted_labels"], adata.obs.columns)
        adata.obs.loc[:, cols].to_csv(tsv, sep="\t")
    sc.pl.umap(adata, color=leiden, save="_%s_beside.png" % leiden)
    sc.pl.umap(adata, color=leiden, save="_%s_ondata.png" % leiden, legend_loc="on data", legend_fontsize=2)
    if dotplot is not None:
        sc.pl.dotplot(adata, var_names=dotplot, groupby=leiden, save="%s.png" % leiden, standard_scale="var")
    for vv in np.intersect1d(["pct_counts_mt", "doublet_score", "log1p_total_counts"], adata.obs.columns):
        sc.pl.violin(adata, vv, groupby=leiden, save="_%s_%s.png" % (leiden, vv))
    if use_rgg:
        sc.tl.dendrogram(adata, groupby=leiden)
        with sw("Ranking genes"):
            sc.tl.rank_genes_groups(adata, groupby=leiden, method="wilcoxon", pts=True)
        sc.pl.rank_genes_groups_dotplot(adata, save="rgg_%s.png" % leiden, n_genes=rgg_ng)
        sc.pl.rank_genes_groups_matrixplot(adata, save="rgg_%s.png" % leiden, n_genes=rgg_ng)
        sc.pl.rank_genes_groups_heatmap(adata, save="_rgg_%s.png" % leiden, n_genes=rgg_ng)
    if output is not None:
        with sw("Re-setting counts") as _:
            adata.X = adata.layers["raw"].copy()
            del adata.layers["raw"]
        with sw("Writing to H5AD"):
            adata.write_h5ad(output, compression="gzip", compression_opts=compression)
    return adata

if __name__ == "__main__":
    import benj
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", dest="h5ad", required=True)
    ap.add_argument("-o", "--output")
    ap.add_argument("-t", "--tsv")
    ap.add_argument("-b", "--batch", type=str, default=None)
    ap.add_argument("-p", "--plot", nargs="+")
    ap.add_argument("-r", "--resolution", default=1., type=float)
    ap.add_argument("-l", "--leiden", default="overall_clust")
    ap.add_argument("--prefix", default="C")
    ap.add_argument("--hvg", default=0, type=int)
    ap.add_argument("--no-use-combat", dest="use_combat", action="store_false")
    ap.add_argument("--use-combat", dest="use_combat", action="store_true")
    ap.add_argument("--no-use-harmony", dest="use_harmony", action="store_false")
    ap.add_argument("--use-harmony", dest="use_harmony", action="store_true")
    ap.add_argument("--no-use-bbknn", dest="use_bbknn", action="store_false")
    ap.add_argument("--use-bbknn", dest="use_bbknn", action="store_true")
    ap.add_argument("--no-use-scaling", dest="use_scaling", action="store_false")
    ap.add_argument("--use-scaling", dest="use_scaling", action="store_true")
    ap.add_argument("--dotplot", nargs="+")
    ap.add_argument("--target-sum", type=int, default=0)
    ap.add_argument("--celltypist")
    ap.add_argument("--compression", type=int, default=9)
    ap.add_argument("--min-dist", type=float, default=0.5)
    ap.add_argument("--no-rank-genes", dest="use_rgg", action="store_false")
    ap.add_argument("--rank-genes", dest="use_rgg", action="store_true")
    ap.add_argument("--max-iter-harmony", type=int, default=50)
    ap.set_defaults(use_combat=False, use_harmony=False, use_bbknn=True, use_scaling=False, use_rgg=True)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    adata = benj.parse_anndata(**args)
    integrate(adata, **args)
