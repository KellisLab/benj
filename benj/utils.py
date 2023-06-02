
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
