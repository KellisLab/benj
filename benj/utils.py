
def convert_X(X, dtype_tbl=None):
    import numpy as np
    import scipy.sparse
    if dtype_tbl is None:
        dtype_tbl = {np.iinfo(np.int8).max: np.int8,
                     np.iinfo(np.uint8).max: np.uint8,
                     np.iinfo(np.int16).max: np.int16,
                     np.iinfo(np.uint16).max: np.uint16,
                     np.iinfo(np.int32).max: np.int32,
                     np.iinfo(np.uint32).max: np.uint32,
                     np.iinfo(np.int64).max: np.int64}
    if scipy.sparse.issparse(X):
        D = X.data
    else:
        D = X
    if isinstance(X.dtype, (int, np.integer)):
        data_max = np.abs(D).max()
        for dtype_max in sorted(dtype_tbl.keys()):
            if data_max < dtype_max:
                return X.astype(dtype_tbl[dtype_max])
    elif np.allclose(D, np.round(D)):
        data_max = np.abs(D).max()
        for dtype_max in sorted(dtype_tbl.keys()):
            if data_max < dtype_max:
                return X.astype(dtype_tbl[dtype_max])
    return X

def is_norm_log(adata, target_sum=10000):
    from scipy.special import logsumexp
    import numpy as np
    if adata.dtype.kind in ["u", "i"]:
        return False
    X = adata.chunk_X()
    lse = logsumexp(X, axis=0)
    total = np.exp(lse) - X.shape[1]
    return np.allclose(total, target_sum)

def index_of(values, idx):
    import numpy as np
    sort_idx = np.argsort(idx)
    values_idx = sort_idx[np.searchsorted(idx, values, sorter=sort_idx)]
    return values_idx

def weighted_pearson_correlation(A, B, wt=None):
    import numpy as np
    if wt is None:
        wt = np.ones(A.shape[1])
    wt = np.ravel(wt) / np.sum(wt)
    Am = (A.T - np.dot(A, wt)).T
    Bm = (B.T - np.dot(B, wt)).T
    Am_sum_sq = np.dot(Am * Am, wt)
    Bm_sum_sq = np.dot(Bm * Bm, wt)
    numer = np.dot(Am * Bm, wt)
    denom = np.sqrt(Am_sum_sq * Bm_sum_sq)
    cor = np.divide(numer, denom, out=np.zeros_like(numer), where=denom != 0)
    return cor

def pseudobulk_valid_columns(obs, inv):
    ### find cols that don't change within cols
    from functools import reduce
    import numpy as np
    import pandas as pd
    from tqdm.auto import tqdm
    goodcols = set(obs.columns.values)
    for i in tqdm(np.arange(np.max(inv)+1)):
        df = obs.loc[i == inv, list(goodcols)]
        flag = df.apply(pd.Series.nunique, axis=0) == 1
        goodcols &= set(df.columns.values[flag])
    return list(goodcols)

def pseudobulk(adata, cols, which:str="sum", dense:bool=True):
    import numpy as np
    import scipy.sparse
    import pandas as pd
    from functools import reduce
    from tqdm.auto import tqdm
    import anndata
    ug, gidx, ginv, gcnt = np.unique(adata.obs.groupby(cols).ngroup(), return_inverse=True, return_counts=True, return_index=True)
    cols = pseudobulk_valid_columns(adata.obs, ginv)
    if which == "sum":
        S = scipy.sparse.csr_matrix((np.ones_like(ginv), (ginv, np.arange(len(ginv)))),
                                    shape=(len(ug), adata.shape[0]))
    elif which == "mean":
        S = scipy.sparse.csr_matrix(((1/gcnt)[ginv], (ginv, np.arange(len(ginv)))),
                                    shape=(len(ug), adata.shape[0]))
    else:
        raise ValueError("Which is not supported")
    obs = adata.obs.loc[:, cols].iloc[gidx, :]
    if adata.isbacked:
        if dense:
            PX = np.zeros(shape=(S.shape[0], adata.shape[1]), dtype=S.dtype)
        else:
            PX = scipy.sparse.lil_matrix((S.shape[0], adata.shape[1]), dtype=S.dtype)
        for X, start, end in tqdm(adata.chunked_X()):
            X = scipy.sparse.coo_matrix(X)
            X = scipy.sparse.csr_matrix((X.data, (X.row + start, X.col)), adata.shape)
            PX += S.dot(X)
        if dense:
            PX = np.asarray(PX)
        elif scipy.sparse.issparse(PX):
            PX = PX.tocsr()
    else:
        PX = S.dot(adata.X)
        if dense and scipy.sparse.issparse(PX):
            PX = np.asarray(PX.todense())
    return anndata.AnnData(X=convert_X(PX),
                           obs=obs,
                           var=adata.var,
                           layers={k: convert_X(S.dot(adata.layers[k])) for k in adata.layers.keys()})
