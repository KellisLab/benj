
def _extend_pre_ranges(df, upstream:int=0, downstream:int=0, start:str="Start", end:str="End", strand:str="Strand"):
    df.loc[df[strand] == "+", start] -= upstream
    df.loc[df[strand] == "-", start] -= downstream
    df.loc[df[strand] == "+", end] += downstream
    df.loc[df[strand] == "-", end] += upstream
    return df

def estimate_genes_archr(adata, gtf:str,
                         min_upstream:int=1000,
                         max_upstream:int=100000,
                         min_downstream:int=1000,
                         max_downstream:int=100000,
                         gene_upstream:int=5000,
                         gene_downstream:int=0,
                         target_sum:int=1000,
                         gene_scale_factor:float=5.,
                         peak_column:str=None, ## If not provided, will use peak index
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
    with sw("Extending ranges"):
        gf["gene_length"] = gf["End"] - gf["Start"]
        ### scale by inverse gene length
        gs = 1/gf["gene_length"].values
        ### min max scale
        gs = (gs - np.min(gs)) / (np.max(gs) - np.min(gs))
        ### expand to 1 to GSF range. But take log, so log1p((gs-1)*s) = log(1 + (gs - 1) * s)
        gf["log_gene_scale"] = np.log1p((gene_scale_factor - 1) * gs)
        gf.index = gf["gene_id"].values
        gf = _extend_pre_ranges(gf, upstream=gene_upstream, downstream=gene_downstream)
        gf["MinStart"] = gf["Start"].values
        gf["MinEnd"] = gf["End"].values
        gf = _extend_pre_ranges(gf, upstream=min_upstream, downstream=min_downstream, start="MinStart", end="MinEnd")
        gr = pyranges.from_dict({"Chromosome": gf["Chromosome"],
                                 "Start": gf["Start"],
                                 "End": gf["End"],
                                 "gene_id": gf["gene_id"]})
        mingr = pyranges.from_dict({"Chromosome": gf["Chromosome"],
                                    "Start": gf["MinStart"],
                                    "End": gf["MinEnd"],
                                    "gene_id": gf["gene_id"]})
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
        ## Once peak ranges are gathered, find intersecting gene bodies:
        inter_df = pr.join(gr).df.loc[:, ["peak_name", "gene_id"]]
        inter_df["Distance"] = 0
        ## Then, find genes with minimum upstream/downstream distance away
        min_df = pr.join(mingr).df.loc[:, ["peak_name", "gene_id"]]
        ## diff. is accurate, unless overlapping intervals. Do not need to worry, as duplicates from inter_df will take care
        diff = pf.loc[min_df["peak_name"].values, ["Start", "Start", "End", "End"]].values.astype(int) - gf.loc[min_df["gene_id"].values, ["Start", "End", "Start", "End"]].values.astype(int)
        min_df["Distance"] = np.abs(diff).min(1) + 1
        ## Finally, find distal. Only need nearest gene
        distance_df = pr.nearest(gr).df.loc[:, ["peak_name", "gene_id", "Distance"]]
        ## Concat such that 1) prioritized intersections, then 2) minimum distance away, then 3) distal
        df = pd.concat([inter_df, min_df, distance_df]).drop_duplicates(["peak_name", "gene_id"])
        df["weight"] = np.exp(-1 - np.abs(df["Distance"]) / 5000. + gf.loc[df["gene_id"].values, "log_gene_scale"].values)
    with sw("Calculating accessibility"):
        S = scipy.sparse.csr_matrix((df["weight"].values,
                                     (pf.index.get_indexer(df["peak_name"].values),
                                      gf.index.get_indexer(df["gene_id"].values))),
                                    shape=(pf.shape[0], gf.shape[0]))
        if layer is not None and layer in adata.layers:
            X = adata.layers[layer]
        else:
            X = adata.X
        gdata = anndata.AnnData(X.dot(S), obs=adata.obs, var=gf.loc[:, ["gene_id", "gene_name"]], dtype=np.float32, obsm=adata.obsm)
    gdata.var_names = gdata.var["gene_name"].values
    gdata.var_names_make_unique()
    del gdata.var["gene_name"]
    gdata.var.columns = ["gene_ids"]
    sc.pp.normalize_total(gdata, target_sum=target_sum)
    return gdata
