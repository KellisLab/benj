
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
