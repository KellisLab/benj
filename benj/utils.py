
from typing import Union, List, Dict
def pg(adata, color_map="Reds"):
    import scanpy as sc
    return lambda x, **kwargs: sc.pl.umap(adata, color=x, color_map=color_map, save="_%s.png" % x, **kwargs)


def read_elems(path: str, elems: Union[str, List[str]]) -> Union[Dict[str, any], any]:
    import h5py
    import anndata.experimental
    with h5py.File(path, "r") as F:
        if isinstance(elems, str):
            return anndata.experimental.read_elem(F[elems])
        else:
            return {elem: anndata.experimental.read_elem(F[elem]) for elem in elems}

def filter_LSI(adata, qc_cols, cor_cutoff:float=0.8, sw=None):
    import numpy as np
    if sw is None:
        from .template import stopwatch
        sw = stopwatch()
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

def convert_X(X, dtype_tbl=None):
    import numpy as np
    import scipy.sparse
    import gc
    if dtype_tbl is None:
        dtype_tbl = {np.iinfo(np.int8).max: np.int8,
                     np.iinfo(np.uint8).max: np.uint8,
                     np.iinfo(np.int16).max: np.int16,
                     np.iinfo(np.uint16).max: np.uint16,
                     np.iinfo(np.int32).max: np.int32,
                     np.iinfo(np.uint32).max: np.uint32,
                     np.iinfo(np.int64).max: np.int64}
    if scipy.sparse.issparse(X):
        X.sum_duplicates()
        D = X.data
    else:
        D = X
    if isinstance(X.dtype, (int, np.integer)):
        if len(D) == 0:
            data_max = 0
        else:
            data_max = np.abs(D).max()
        for dtype_max in sorted(dtype_tbl.keys()):
            if data_max <= dtype_max:
                return X.astype(dtype_tbl[dtype_max])
    elif np.allclose(D, np.round(D)):
        if len(D) == 0:
            data_max = 0
        else:
            data_max = np.abs(D).max()
        for dtype_max in sorted(dtype_tbl.keys()):
            if data_max <= dtype_max:
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

def anndata_batched_generator(adata, batch_size:int=10000):
    from anndata import AnnData
    from anndata.experimental import AnnCollection, AnnCollectionView
    if isinstance(adata, AnnData):
        adata = AnnCollection({"data": adata}, join_vars="inner")
    return adata.iterate_axis(batch_size)

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

def dotplot_summarize(adata, groupby, expression_cutoff:float=0, which="mean"):
    import scanpy as sc
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    adata.layers["expressed"] = adata.X > expression_cutoff
    return pseudobulk(adata, cols=groupby, which=which)
