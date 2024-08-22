
from typing import Union, List, Dict
def pg(adata, color_map="Reds"):
    import scanpy as sc
    return lambda x, **kwargs: sc.pl.umap(adata, color=x, color_map=color_map, save="_%s.png" % x, **kwargs)

def read_elems(path: str, elems: Union[str, List[str]], retry:int=2) -> Union[Dict[str, any], any]:
    import h5py
    import anndata.experimental
    if retry < 0:
        return {}
    try:
        with h5py.File(path, "r") as F:
            if isinstance(elems, str):
                return anndata.experimental.read_elem(F[elems])
            else:
                return {elem: anndata.experimental.read_elem(F[elem]) for elem in elems}
    except KeyError:
        import time
        print("Key error for %s, sleeping 5 seconds and retrying..." % path)
        time.sleep(5)
        ### NFS error
        return read_elems(path=path, elems=elems, retry=retry-1)

def filter_LSI(adata, qc_cols, cor_cutoff:float=0.8, sw=None):
    import numpy as np
    if sw is None:
        from .timer import template as stopwatch
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

def convert_X(X, dtype_tbl=None, signed_only:bool=True):
    import numpy as np
    import scipy.sparse
    import gc
    if dtype_tbl is None:
        if signed_only:
            dtype_tbl = {np.iinfo(np.int8).max: np.int8,
                         np.iinfo(np.int16).max: np.int16,
                         np.iinfo(np.int32).max: np.int32,
                         np.iinfo(np.int64).max: np.int64}
        else:
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

def as_ranges(var, interval="interval", extend_left:int=0, extend_right:int=0, index_col="_index"):
    import numpy as np
    import pyranges
    rvar = var[interval].str.extract(r"^([^:]+):([0-9]+)-([0-9]+)$")
    var = var.loc[rvar.isna().sum(1) == 0, :].copy()
    rvar = rvar.loc[rvar.isna().sum(1) == 0, :]
    var["Chromosome"] = rvar[0]
    var["Start"] = rvar[1].astype(int) - extend_left
    var["End"] = rvar[2].astype(int) + extend_right
    var[index_col] = var.index.values
    return pyranges.from_dict({k: var[k].values for k in var.columns})

def powerlaw_at_distance(distance, hic_power=0.87):
    """ABC raw distance; has to be scaled by denom"""
    import numpy as np
    scale = -4.80 + 11.63 * hic_power
    offset = np.clip(np.abs(distance), 5000, np.Inf)
    return np.exp(scale + -1 * hic_power * np.log1p(offset))

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

def is_nonzero_combinatorial(adata, genes, label:str="Combination", prefix:str="Combo ", inplace:bool=True):
    import numpy as np
    import pandas as pd
    import scipy.sparse
    I = adata.var_names.get_indexer(genes)
    X = adata.X[:, I]
    if scipy.sparse.issparse(X):
        X = X.todense().A
    B = (X > X.mean(0)) @ [1 << x for x in range(X.shape[1])]
    combos = []
    for i in range(len(set(B))):
        subset = []
        for j, g in enumerate(genes):
            if i & (1 << j) == (1 << j):
                subset.append(g)
        if subset:
            combos.append(" ".join(subset))
        else:
            combos.append("None")
    df = pd.DataFrame({label: np.asarray(combos)[B]}, index=adata.obs_names)
    xf = pd.get_dummies(df[label]).map(lambda x: 1 if x else np.nan)
    xf.columns = ["%s%s" % (prefix, cn) for cn in xf.columns]
    df = pd.concat((df, xf), axis=1)
    if inplace:
        for cn in df.columns:
            adata.obs[cn] = df[cn]
    else:
        return df

    
def is_gzip(filename):
    with open(filename, "rb") as f:
        return f.read(2) == b"\x1f\x8b"

def parse_gmt(filename, genes):
    import gzip
    import pandas as pd
    import scipy.sparse
    import anndata
    genes = pd.Index(genes)
    mat_rows = []
    with (gzip.open(filename, "r") if is_gzip(filename) else open(filename, "r")) as F:
        for row_idx, line in enumerate(F):
            tokens = line.decode("utf-8").strip().split("\t")
            gsname = tokens[0]
            ind = genes.get_indexer(tokens[2:])
            val_ind = ind[ind != -1]
            row = scipy.sparse.lil_matrix((1, len(genes)), dtype="i1")
            row[0, val_ind] = 1
            mat_rows.append(row)
    mat = scipy.sparse.lil_matrix((len(mat_rows), len(genes)), dtype="i1")
    for i, row in enumerate(mat_rows):
        mat[i, :] = row
    return anndata.AnnData(mat.T.tocsr(), obs=pd.DataFrame(index=genes))

def dotplot_summarize(adata, groupby, expression_cutoff:float=0, which="mean"):
    import scanpy as sc
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    adata.layers["expressed"] = adata.X > expression_cutoff
    return pseudobulk(adata, cols=groupby, which=which)

def ac_var_qc(ac, batch_size:int=10000, meansd:bool=True, layer:str=None):
    import numpy as np
    import pandas as pd
    import anndata
    import scanpy as sc
    total_counts = np.zeros(ac.shape[1], dtype="i8")
    n_cells = 0
    mean = np.zeros(ac.shape[1], dtype="f8")
    var = np.zeros(ac.shape[1], dtype="f8")
    for batch, _ in tqdm(ac.iterate_axis(10000, shuffle=meansd)):
        if layer is not None:
            X = batch.X
        else:
            X = batch.layers[layer]
        total_counts += np.ravel(X.sum(0))
        if meansd:
            batch = anndata.AnnData(X)
            sc.pp.normalize_total(batch, target_sum=10000)
            sc.pp.log1p(batch)
            s1 = np.ravel(batch.X.sum(0))
            s2 = np.ravel(batch.X.multiply(batch.X).sum(0))
            mean_new = (n_cells * mean + s1) / (n_cells + X.shape[0])
            ### Young and Cramer update
            var_new = (n_cells * (var + mean*mean) + s2) / (n_cells + X.shape[0]) - mean_new*mean_new
            mean = mean_new
            var = var_new
        n_cells += X.shape[0]
    var = pd.DataFrame(index=ac.var_names)
    var["total_counts"] = total_counts
    if meansd:
        var["mean"] = mean
        var["std"] = np.sqrt(var)
    return var

def transform_matrix_index(index_from, index_to):
    import scipy.sparse
    import numpy as np
    import pandas as pd
    if not isinstance(index_from, pd.Index):
        index_from = pd.Index(index_from)
    if not isinstance(index_to, pd.Index):
        index_to = pd.Index(index_to)
    comm = index_from.intersection(index_to)
    S = scipy.sparse.csr_matrix((np.ones(len(comm), dtype=int),
                                 (index_from.get_indexer(comm),
                                  index_to.get_indexer(comm))),
                                 shape=(len(index_from), len(index_to)), dtype=int)
    return S

def leiden_multiplex(mdata, resolution:float=1., key_added="mleiden",
                     neighbors_key:str=None, n_iterations=-1,
                     prefix="C", **kwargs):
    import numpy as np
    import scipy.sparse
    from scanpy._utils import get_igraph_from_adjacency
    from scanpy.tools._utils import _choose_graph
    import leidenalg
    obs_names = mdata.obs_names
    tbl = {}
    for k in mdata.mod.keys():
        ### Translate per-modality indices into union of anndata indices
        S = scipy.sparse.csr_matrix((np.ones(mdata.mod[k].shape[0], dtype=int),
                                     (obs_names.get_indexer(mdata.mod[k].obs_names),
                                      np.arange(mdata.mod[k].shape[0]))),
                                     shape=(len(obs_names), mdata.mod[k].shape[0]),
                                    dtype=int)
        ### Sandwich to become UxU matrix
        adjacency = _choose_graph(mdata.mod[k], None, neighbors_key)
        g = S.dot(adjacency.dot(S.T))
        tbl[k] = get_igraph_from_adjacency(g)
    clust, _ = leidenalg.find_partition_multiplex(tbl.values(),
                                                  leidenalg.RBConfigurationVertexPartition,
                                                  resolution_parameter=resolution,
                                                  n_iterations=n_iterations,
                                                  **kwargs)
    mdata.obs[key_added] = ["%s%s" % (prefix, s) for s in np.asarray(clust).astype(str)]
