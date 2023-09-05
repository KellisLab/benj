
def _extend_pre_ranges(df, upstream:int=0, downstream:int=0, start:str="Start", end:str="End", strand:str="Strand"):
    strand_rm = False
    if strand not in df.columns:
        strand_rm = True
        df[strand] = "+"
    df.loc[df[strand] == "+", start] -= upstream
    df.loc[df[strand] == "-", start] -= downstream
    df.loc[df[strand] == "+", end] += downstream
    df.loc[df[strand] == "-", end] += upstream
    if strand_rm:
        del df[strand]
    return df

def _extend_ranges(gf,
                   min_upstream:int=1000,
                   max_upstream:int=100000,
                   min_downstream:int=1000,
                   max_downstream:int=100000,
                   gene_upstream:int=5000,
                   gene_downstream:int=0,
                   gene_scale_factor:float=5.):
    import numpy as np
    gf["gene_length"] = gf["End"] - gf["Start"]
    ### In case of duplicated indices, sum lengths across indices
    gb = gf.groupby(level=0).agg(gene_length=("gene_length", sum))
    gf["gene_length"] = gb.loc[gf.index.values, "gene_length"].values
    ### scale by inverse gene length
    gs = 1/gf["gene_length"].values
    ### min max scale, plus epsilon to avoid divide by 0
    gs = (gs - np.min(gs)) / (np.max(gs) - np.min(gs)) + 1e-300
    ### expand to 1 to GSF range. But take log, so log1p((gs-1)*s) = log(1 + (gs - 1) * s)
    gf["log_gene_scale"] = np.log1p((gene_scale_factor - 1) * gs)
    #gf.index = gf["gene_id"].values
    gf = _extend_pre_ranges(gf, upstream=gene_upstream, downstream=gene_downstream)
    gf["MinStart"] = gf["Start"].values
    gf["MinEnd"] = gf["End"].values
    gf = _extend_pre_ranges(gf, upstream=min_upstream, downstream=min_downstream, start="MinStart", end="MinEnd")
    gf["interval"] = ["%s:%d-%d" % (chrom, start, end) for chrom, start, end in zip(gf["Chromosome"], gf["Start"], gf["End"])]
    gf["min_interval"] = ["%s:%d-%d" % (chrom, start, end) for chrom, start, end in zip(gf["Chromosome"], gf["MinStart"], gf["MinEnd"])]
    return gf

def estimate_features_archr(adata, feature_df,
                            min_upstream:int=1000,
                            max_upstream:int=100000,
                            min_downstream:int=1000,
                            max_downstream:int=100000,
                            gene_upstream:int=5000,
                            gene_downstream:int=0,
                            target_sum:int=None,
                            gene_scale_factor:float=5.,
                            peak_column:str=None, ## If not provided, will use peak index
                            feature_column:str=None,
                            var_column_tolerance:float=0.999, ### Tolerance used for determining if a feature should be kept, in case of weird things
                            distal:bool=True, ### Use nearest gene to a peak if unassigned
                            log1p:bool=False,
                            layer:str=None):
    import numpy as np
    import pandas as pd
    import scipy.sparse
    import pyranges
    import anndata
    import scanpy as sc
    from .timer import template
    sw = template()
    if not isinstance(feature_df, pd.DataFrame):
        raise ValueError("Feature_df is not a dataframe")
    if not np.all(np.isin(["Chromosome", "Start", "End"], feature_df.columns)):
        raise ValueError("Feature_df does not have ranges")
    if feature_column is not None:
        feature_df.index = feature_df[feature_column].values
    with sw("Extending ranges"):
        gf = _extend_ranges(feature_df, min_upstream=min_upstream, min_downstream=min_downstream,
                            max_upstream=max_upstream, max_downstream=max_downstream,
                            gene_upstream=gene_upstream, gene_downstream=gene_downstream,
                            gene_scale_factor=gene_scale_factor)
        gr = pyranges.from_dict({"Chromosome": gf["Chromosome"],
                                 "Start": gf["Start"],
                                 "End": gf["End"],
                                 "feature_index": gf.index.values,
                                 "feature_interval": gf["interval"].values})
        mingr = pyranges.from_dict({"Chromosome": gf["Chromosome"],
                                    "Start": gf["MinStart"],
                                    "End": gf["MinEnd"],
                                    "feature_index": gf.index.values,
                                    "feature_interval": gf["interval"].values})
    ##
    ## Now, get peak ranges (pr) from peak frame (pf)
    ##
    with sw("Extracting peak ranges"):
        if peak_column is None:
            pstr = adata.var_names.values
        else:
            pstr = adata.var[peak_column].values
        pf = pd.DataFrame([x.replace(":", "-", 1).split("-") for x in pstr], columns=["Chromosome", "Start", "End"], index=adata.var_names)
        pr = pyranges.from_dict({"Chromosome": pf["Chromosome"], "Start": pf["Start"], "End": pf["End"], "peak_name": pf.index.values})
    with sw("Calculating overlaps"):
        iif = gf.drop_duplicates("interval")
        iif.index = iif["interval"]
        ## Once peak ranges are gathered, find intersecting gene bodies:
        inter_df = pr.join(gr).df.loc[:, ["peak_name", "feature_index", "feature_interval"]]
        inter_df["Distance"] = 0
        ## Then, find genes with minimum upstream/downstream distance away
        min_df = pr.join(mingr).df.loc[:, ["peak_name", "feature_index", "feature_interval"]]
        ## diff. is accurate, unless overlapping intervals. Do not need to worry, as duplicates from inter_df will take care
        diff = pf.loc[min_df["peak_name"].values, ["Start", "Start", "End", "End"]].values.astype(int) - iif.loc[min_df["feature_interval"].values, ["Start", "End", "Start", "End"]].values.astype(int)
        min_df["Distance"] = np.abs(diff).min(1) + 1
        ## Finally, find distal. Only need nearest gene

        if distal:
            distance_df = pr.nearest(gr).df.loc[:, ["peak_name", "feature_index", "feature_interval", "Distance"]]
            ## Concat such that 1) prioritized intersections, then 2) minimum distance away, then 3) distal
            df = pd.concat([inter_df, min_df, distance_df]).drop_duplicates(["peak_name", "feature_index"])
        else:
            df = pd.concat([inter_df, min_df]).drop_duplicates(["peak_name", "feature_index"])
        df["weight"] = np.exp(-1 - np.abs(df["Distance"]) / 5000. + iif.loc[df["feature_interval"].values, "log_gene_scale"].values)
    with sw("Calculating accessibility"):
        if gf.index.duplicated().sum() > 0:
            ### Get columns that are the same across repeated indices
            nf = gf.groupby(level=0).nunique()
            gf = gf.loc[~gf.index.duplicated(keep="first"), nf.columns[(nf==1).mean() >= var_column_tolerance]]
        S = scipy.sparse.csr_matrix((df["weight"].values,
                                     (pf.index.get_indexer(df["peak_name"].values),
                                      gf.index.get_indexer(df["feature_index"].values))),
                                    shape=(pf.shape[0], gf.shape[0]))
        if layer is not None and layer in adata.layers:
            X = adata.layers[layer]
        else:
            X = adata.X
        gdata = anndata.AnnData(X.dot(S), obs=adata.obs, var=gf, dtype=np.float32, obsm=adata.obsm, obsp=adata.obsp,
                                uns={k: v for k, v in adata.uns.items() if k in ["neighbors", "files", "lsi", "pca", "umap", "leiden"]})
    if target_sum is not None and target_sum > 0:
        sc.pp.normalize_total(gdata, target_sum=target_sum)
    else:
        print("Using median normalization")
        sc.pp.normalize_total(gdata)
    if log1p:
        sc.pp.log1p(gdata)
    return gdata

def estimate_genes_archr(adata, gtf:str,
                         min_upstream:int=1000,
                         max_upstream:int=100000,
                         min_downstream:int=1000,
                         max_downstream:int=100000,
                         gene_upstream:int=5000,
                         gene_downstream:int=0,
                         target_sum:int=None,
                         gene_scale_factor:float=5.,
                         peak_column:str=None, ## If not provided, will use peak index
                         feature_column:str="gene_id", ### If not provided, will use feature index
                         log1p:bool=False,
                         layer:str=None):
    import numpy as np
    import pandas as pd
    import scipy.sparse
    import pyranges
    import anndata
    import scanpy as sc
    from .timer import template
    sw = template()
    with sw("Reading GTF"):
        gf = pyranges.read_gtf(gtf).df
        gf = gf.loc[gf["Feature"] == "gene", :]
        if feature_column in gf.columns:
            gf.index = gf[feature_column].values.astype(str)
    gdata = estimate_features_archr(adata, feature_df=gf,
                                   min_upstream=min_upstream, min_downstream=min_downstream,
                                   max_upstream=max_upstream, max_downstream=max_downstream,
                                   gene_upstream=gene_upstream, gene_downstream=gene_downstream,
                                   target_sum=target_sum, gene_scale_factor=gene_scale_factor,
                                   peak_column=peak_column,
                                   feature_column=feature_column,
                                   log1p=log1p, layer=layer)
    gdata.var = gdata.var.loc[:, ["gene_id", "gene_name"]]
    gdata.var_names = gdata.var["gene_name"].values
    gdata.var_names_make_unique()
    del gdata.var["gene_name"]
    gdata.var.columns = ["gene_ids"]
    add_gene_length(gdata.var, gtf)
    return gdata

def get_tss(tss:str):
    import pandas as pd
    tss = pd.read_csv(tss, sep="\t", header=None)
    tss.columns = ["Chromosome", "Start", "End", "gene_id", "score", "strand"]
    df = tss.groupby(["Chromosome", "gene_id", "strand"]).agg(left=("Start", "min"),
                                                              right=("End", "max")).reset_index()
    df["interval"] = df["Chromosome"] + ":" + df["left"].astype(str) + "-" + df["right"].astype(str)
    df.index = df["gene_id"].values
    return df

def add_interval(var, tss:str, inplace=True):
    tf = get_tss(tss)
    interval = [tf["interval"].get(g, "NA") for g in var["gene_ids"]]
    if inplace:
        var["interval"] = interval
    else:
        import pandas as pd
        return pd.Series(interval, index=var.index, name="interval")

def add_gene_length(var, gtf:str=None, inplace=True):
    import pandas as pd
    import pyranges
    from .timer import template
    sw = template()
    with sw("Reading GTF"):
        gf = pyranges.read_gtf(gtf).df
        gf = gf.loc[gf["Feature"] == "gene",:]
        gf["gene_length"] = gf["End"] - gf["Start"]
        gf.index = gf["gene_id"].values
    gl = [gf["gene_length"].get(g, -1) for g in var["gene_ids"]]
    gs = [gf["Strand"].get(g, "*") for g in var["gene_ids"]]
    if not inplace:
        var = var.copy()
    var["gene_length"] = gl
    var["strand"] = gs
    if not inplace:
        return var.loc[:, ["gene_length", "strand"]]

def add_gene_info(var, gene_info:str=None, inplace=True):
    """Add a STAR geneInfo.tab file to .var"""
    import pandas as pd
    df = pd.read_csv(gene_info, sep="\t", skiprows=[0], header=None)
    df.index = df[0].values
    gi = [df[2].get(g, "NA") for g in var["gene_ids"]]
    if inplace:
        var["gene_type"] = gi
    else:
        return pd.Series(gi, index=var.index, name="gene_type")
