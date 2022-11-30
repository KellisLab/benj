
def idf(adata, features, chunk_size):
    import numpy as np
    n, m = adata.shape
    count = np.zeros(m, dtype=np.uint64)
    for batch, _, _ in adata.chunked_X(chunk_size):
        uidx, idx_count = np.unique(batch.nonzero()[1], return_counts=True)
        count[uidx] = count[uidx] + idx_count
    if features is not None:
        count = count[features]
    return np.log1p(n / (1 + count))


def _chunked_X(adata, select, replace=True):
    """Use anndata.chunked_X, but instead of densifying,
    keep matrix as-is."""
    import numpy as np
    reverse = None
    select = select if select < adata.n_obs else adata.n_obs
    choice = np.random.choice(adata.n_obs, select, replace)
    if adata.isbacked:
        indices, reverse = np.unique(choice, return_inverse=True)
        selection = adata.X[indices.tolist()]
    else:
        selection = adata.X[choice]
    return selection if reverse is None else selection[reverse]

def spectral(adata, n_comps=50, features="selected",
             random_state=0,
             sample_size=0.25,
             chunk_size=20000,
             distance_metric="jaccard",
             inplace=True):
    """Use SnapATAC2 spectral, but without the rust interface"""
    from snapatac2.tl._spectral import Spectral
    import numpy as np
    from tqdm import tqdm
    import scipy.sparse
    np.random.seed(random_state)
    if features in adata.var.columns:
        features = adata.var[features].values
    else:
        print("Warning: using all features")
        features = np.repeat(True, adata.shape[1])
    weights = idf(adata, features=features, chunk_size=chunk_size)
    n_comps = min(adata.shape[0], adata.shape[1], n_comps)
    sample_size = int(adata.shape[0] * sample_size)
    ### Get sample
    S = _chunked_X(adata, sample_size, replace=True)
    if not scipy.sparse.issparse(S):
        S = scipy.sparse.csr_matrix(S)
    if features is not None:
        S = S[:, features]
    model = Spectral(n_comps, distance_metric, weights)
    model.fit(S.astype(np.float64))
    del S
    print("Perform Nystroem extension")
    for batch, _, _ in adata.chunked_X(chunk_size):
        if distance_metric == "jaccard":
            batch = (batch > 0)
        if features is not None:
            batch = batch[:, features]
        model.extend(batch.astype(np.float64))
    result = model.transform()
    if inplace:
        adata.uns['spectral_eigenvalue'] = result[0]
        adata.obsm['X_spectral'] = np.real(result[1])
    else:
        return result[0], np.real(result[1])


def n_spectral(eigenvalues, S=1.0, curve="convex", direction="decreasing"):
    import numpy as np
    from kneed import KneeLocator
    kl = KneeLocator(np.arange(len(eigenvalues)),
                     eigenvalues, S=S, direction=direction,
                curve=curve)
    return kl.knee

def find_features(adata, resolution=1., features="selected", n_features=10000, blacklist="blacklist", min_comps=5, chunk_size=20000, max_val=4):
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import anndata
    import scipy.sparse
    from tqdm.auto import tqdm
    n_pcs = max(n_spectral(adata.uns["spectral_eigenvalue"]), min_comps)
    sc.pp.neighbors(adata, n_pcs=n_pcs, use_rep="X_spectral")
    sc.tl.leiden(adata, resolution=resolution)
    ### pseudobulk per cluster
    ul, linv = np.unique(adata.obs["leiden"], return_inverse=True)
    P = np.zeros((len(ul), adata.shape[1]), dtype=int)
    for X, begin, end in adata.chunked_X(chunk_size=chunk_size):
        if scipy.sparse.issparse(X):
            X.data = np.clip(X.data, 0, max_val)
        else:
            X = np.clip(X, 0, max_val)
        this_L = linv[begin:end]
        S = scipy.sparse.csr_matrix((np.ones(len(this_L)), (this_L, np.arange(len(this_L)))),
                                    shape=(len(ul), len(this_L)), dtype=int)
        S = S.dot(X)
        if scipy.sparse.issparse(S):
            S = S.todense()
        P += np.asarray(S)
        del S
    ### then highly variable per cluster
    pdata = anndata.AnnData(P, obs=pd.DataFrame(index=ul), var=adata.var, dtype=int)
    if blacklist is not None and blacklist in pdata.var.columns:
        pdata = pdata[:, ~pdata.var[blacklist]].copy()
    sc.pp.normalize_total(pdata, target_sum=10000)
    sc.pp.log1p(pdata)
    sc.pp.highly_variable_genes(pdata, n_top_genes=n_features, subset=True)
    adata.var[features] = adata.var.index.isin(pdata.var.index.values)

def iterativeSpectral(adata, n_comps=50,
                      features="selected",
                      blacklist="blacklist",
                      n_features=25000,
                      nn_min_comps=5,
                      n_iter=3,
                      resolution=1.,
                      random_state=0,
                      sample_size=0.25,
                      chunk_size=20000,
                      distance_metric="jaccard",
                      inplace=True):
    adata = spectral(adata, n_comps=n_comps, features=features,
                     random_state=random_state,
                     sample_size=sample_size,
                     chunk_size=chunk_size,
                     distance_metric=distance_metric)
    for i in range(1, n_iter):
        print("Iteration", i)
        find_features(adata, resolution=resolution, features=features, blacklist=blacklist, min_comps=nn_min_comps, chunk_size=chunk_size)
        adata = spectral(adata, n_comps=n_comps, features=features,
                         random_state=random_state,
                         sample_size=sample_size,
                         chunk_size=chunk_size,
                         distance_metric=distance_metric)
    return adata
