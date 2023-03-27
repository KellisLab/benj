#!/usr/bin/env python3
def integrate(adata, output=None, batch=None, hvg=0, use_combat=False, use_scaling=False, use_harmony=True, use_bbknn=True, plot=None, leiden="overall_clust", resolution=1., min_dist:float=0.5, dotplot=None, celltypist=None, tsv=None, rgg_ng=5, max_iter_harmony:int=50, prefix="", sw=None, compression:int=9, **kwargs):
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
            adata = adata[adata.obs[batch].isin(cdf.index), :].copy()
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
        sc.pp.normalize_total(adata, target_sum=10000)
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
    if celltypist is not None:
        ct = benj.annotate(adata.raw.to_adata(), majority_voting=True, model=celltypist)
        for cn in ct.predicted_labels.columns:
            adata.obs[cn] = ct.predicted_labels[cn]
        newlabel = "%s_voting" % leiden
        obs = benj.annotate_clusters_from_vote(adata.obs, leiden, "majority_voting", newlabel)
        adata.obs[newlabel] = obs[newlabel]
        del ct, obs
    with sw("Running Leiden"):
        sc.tl.leiden(adata, resolution=resolution, key_added=leiden)
        adata.obs[leiden] = ["%s%s" % (prefix, v) for v in adata.obs[leiden].values.astype(str)]
    if tsv is not None:
        cols = np.intersect1d([leiden, "majority_voting", "predicted_labels"], adata.obs.columns)
        adata.obs.loc[:, cols].to_csv(tsv, sep="\t")
    sc.pl.umap(adata, color=leiden, save="_%s_beside.png" % leiden)
    sc.pl.umap(adata, color=leiden, save="_%s_ondata.png" % leiden, legend_loc="on data")
    if dotplot is not None:
        sc.pl.dotplot(adata, var_names=dotplot, groupby=leiden, save="%s.png" % leiden, standard_scale="var")
    for vv in np.intersect1d(["pct_counts_mt", "doublet_score", "log1p_total_counts"], adata.obs.columns):
        sc.pl.violin(adata, vv, groupby=leiden, save="_%s_%s.png" % (leiden, vv))
    sc.tl.dendrogram(adata, groupby=leiden)
    with sw("Ranking genes"):
        sc.tl.rank_genes_groups(adata, groupby=leiden, method="wilcoxon", pts=True)
    sc.pl.rank_genes_groups_dotplot(adata, save="rgg_%s.png" % leiden, n_genes=rgg_ng)
    sc.pl.rank_genes_groups_matrixplot(adata, save="rgg_%s.png" % leiden, n_genes=rgg_ng)
    sc.pl.rank_genes_groups_heatmap(adata, save="_rgg_%s.png" % leiden, n_genes=rgg_ng)
    with sw("Re-setting counts") as _:
        adata.X = adata.layers["raw"].copy()
    if output is not None:
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
    ap.add_argument("--celltypist")
    ap.add_argument("--compression", type=int, default=9)
    ap.add_argument("--min-dist", type=float, default=0.5)
    ap.add_argument("--max-iter-harmony", type=int, default=50)
    ap.set_defaults(use_combat=False, use_harmony=True, use_bbknn=False, use_scaling=False)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    adata = benj.parse_anndata(**args)
    integrate(adata, **args)
