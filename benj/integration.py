
def integrate_atac(adata, output=None, batch=None, use_harmony:bool=False, use_bbknn:bool=True,
                   leiden="overall_clust", resolution:float=1., prefix="C",
                   tsv=None, min_dist:float=0.3, compression:int=9,
                   min_n_cells_by_counts:int=2, cor_cutoff:float=0.8,
                   max_iter_harmony:int=50,
                   plot=[],
                   qc_cols=["log1p_total_counts"], sw=None, **kwargs):
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import muon as mu
    from muon import atac as ac
    if sw is None:
        from .timer import template as stopwatch
        sw = stopwatch()
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
    if min_n_cells_by_counts > 0 and "n_cells_by_counts" in adata.var.columns:
        with sw("Filtering cells by counts >= %d" % min_n_cells_by_counts):
            mu.pp.filter_var(adata, "n_cells_by_counts", lambda x: x >= min_n_cells_by_counts)
    with sw("Running TF-IDF"):
        adata.X = adata.layers["raw"].copy()
        ac.pp.tfidf(adata)
    with sw("Running LSI"):
        ac.tl.lsi(adata)
        adata.X = adata.layers["raw"].copy()
        del adata.layers["raw"]
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
            sc.external.pp.harmony_integrate(adata, batch, basis=use_rep, adjusted_basis=use_rep_adj, max_iter_harmony=max_iter_harmony)
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
                if ~np.all(adata.obs[col].isna()):
                    sc.pl.umap(adata, color=col, save="_%s.png" % col)
            elif "atac" in adata.uns and "peak_annotation" in adata.uns["atac"] and col in adata.uns["atac"]["peak_annotation"].index:
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
        for col in adata.obs.columns:
            if np.all(adata.obs[col].isna()):
                del adata.obs[col] ### will fail
        adata.write_h5ad(output, compression="gzip", compression_opts=compression)
    return adata


def integrate_rna(adata, output=None, batch=None, hvg:int=0, use_combat:bool=False, use_scaling:bool=False, use_harmony:bool=False, use_bbknn:bool=True, plot=None,
                  leiden:str="overall_clust", resolution:float=1., min_dist:float=0.3,
                  dotplot=None, celltypist=None, tsv=None,
                  rgg_ng:int=5, max_iter_harmony:int=50, prefix:str="C",
                  sw=None, use_rgg:bool=True, target_sum:int=None, compression:int=9, **kwargs):
    import scanpy as sc
    import anndata
    import pandas as pd
    import numpy as np
    if sw is None:
        from .timer import template as stopwatch
        sw = stopwatch()
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
            for col in np.intersect1d(ct.predicted_labels.columns, ["majority_voting", "predicted_labels"]):
                sc.pl.umap(adata, color=col, save="_%s_beside.png" % col)
                sc.pl.umap(adata, color=col, save="_%s_ondata.png" % col, legend_loc="on data", legend_fontsize=2)
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
