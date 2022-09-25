

def annotate(adata, **kwargs):
    import celltypist
    import numpy as np
    import scanpy as sc
    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)
    valid_k = ["over_clustering", "model", "majority_voting"]
    ctargs = {k: v for k, v in kwargs.items() if k in valid_k and v is not None and v != ""}
    ct = celltypist.annotate(adata, **ctargs)
    return ct

def annotate_clusters_from_vote(obs, cluster, annot, newlabel):
    import pandas as pd
    from scipy.stats import hypergeom
    import numpy as np
    cf = obs.groupby([cluster, annot]).count().iloc[:, [0]]
    cf.columns = ["count"]
    cf = cf.loc[cf["count"] > 0, :].reset_index()
    cf = cf.loc[cf[annot] != "", :]
    ct = cf.groupby(cluster).agg(ctotal=("count", sum)).reset_index()
    at = cf.groupby(annot).agg(atotal=("count", sum)).reset_index()
    df = pd.merge(pd.merge(cf, at), ct)
    df["all"] = cf["count"].sum()
    df["mlog10p"] = -hypergeom.logsf(
        df["count"],
        df["all"],
        df["atotal"],
        df["ctotal"]
    ) / np.log(10)
    df = df.sort_values("mlog10p", ascending=False)
    df = df.drop_duplicates(cluster)
    df.index = df[cluster].values
    obs[newlabel] = df.loc[obs[cluster].values, annot].values
    return obs
