
def extract_rank_genes_groups(adata, key="rank_genes_groups_filtered", group_name="group", gene_name="names", pval_cutoff=0.05):
    import pandas as pd
    import numpy as np
    names = pd.DataFrame(adata.uns[key]["names"])
    pvals = pd.DataFrame(adata.uns[key]["pvals_adj"])
    R, C = np.where((~names.isna()) & (pvals <= pval_cutoff))
    df = pd.DataFrame({gene_name: names.values[R, C], group_name: names.columns.values[C]})
    for subkey in adata.uns[key].keys():
        if type(adata.uns[key][subkey]) == np.recarray:
            xf = pd.DataFrame(adata.uns[key][subkey])
            df[subkey] = xf.values[R, C]
    if "pts" in adata.uns[key].keys():
        idx = adata.uns[key]["pts"].index.get_indexer(df[gene_name].values)
        df["pts"] = adata.uns[key]["pts"].values[idx, C]
    return df

def score_genes_by_rgg(adata, rgg):
    import numpy as np
    import scipy.sparse
    import pandas as pd
    import anndata
    df = rgg.loc[adata.var_names.get_indexer(rgg["names"]) >= 0,:].copy()
    df["name_index"] = adata.var_names.get_indexer(df["names"])
    ugrp, grp_index = np.unique(df["group"].values, return_inverse=True)
    df["group_index"] = grp_index
    S = scipy.sparse.csr_matrix((df["scores"].values, (df["name_index"].values, df["group_index"].values)),
                            shape=(adata.shape[1], len(ugrp)))
    #S = S.dot(scipy.sparse.diags(1/np.ravel(S.sum(0))))
    X = S.T.dot(adata.X.T).T
    return anndata.AnnData(X, obs=adata.obs, var=pd.DataFrame(index=ugrp), obsm=adata.obsm, dtype=np.float32)
