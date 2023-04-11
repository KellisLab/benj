
from typing import Union, Optional, Dict
from ._compat import Literal

def extract_rank_genes_groups(adata, key="rank_genes_groups_filtered", group_name="group", gene_name="names", pval_cutoff=0.05):
    import pandas as pd
    import numpy as np
    names = pd.DataFrame(adata.uns[key][gene_name])
    pvals = pd.DataFrame(adata.uns[key]["pvals_adj"])
    R, C = np.where((~names.isna()) & (pvals <= pval_cutoff))
    df = pd.DataFrame({gene_name: names.values[R, C], group_name: names.columns.values[C]})
    for subkey in adata.uns[key].keys():
        if subkey in ["pts", "pts_rest"]:
            idx = adata.uns[key][subkey].index.get_indexer(df[gene_name].values)
            df[subkey] = adata.uns[key][subkey].values[idx, C]
        if isinstance(adata.uns[key][subkey], np.recarray):
            xf = pd.DataFrame(adata.uns[key][subkey])
            df[subkey] = xf.values[R, C]
        elif isinstance(adata.uns[key][subkey], np.ndarray):
            xf = pd.DataFrame(adata.uns[key][subkey])
            df[subkey] = xf.values[R, C]
    return df

def marker_peak_overlap(rgg, mk):
    from sklearn.metrics import jaccard_score
    out = {}
    for k, v in dict(mk).items():
        outk = {}
        for rk, rv in dict(tuple(rgg.groupby("group")["genes"])).items():
            sm = set(v)
            sr = set(rv)
            outk[rk] = len(sm & sr) / len(sm | sr)
        out[k] = outk
    return pd.DataFrame(out)

def score_genes_by_rgg(adata, rgg):
    import numpy as np
    import scipy.sparse
    import pandas as pd
    import anndata
    df = rgg.loc[adata.var_names.get_indexer(rgg["names"]) >= 0, :].copy()
    df["name_index"] = adata.var_names.get_indexer(df["names"])
    ugrp, grp_index = np.unique(df["group"].values, return_inverse=True)
    df["group_index"] = grp_index
    S = scipy.sparse.csr_matrix((df["scores"].values, (df["name_index"].values, df["group_index"].values)),
                            shape=(adata.shape[1], len(ugrp)))
    #S = S.dot(scipy.sparse.diags(1/np.ravel(S.sum(0))))
    X = S.T.dot(adata.X.T).T
    return anndata.AnnData(X, obs=adata.obs, var=pd.DataFrame(index=ugrp), obsm=adata.obsm, dtype=np.float32)
