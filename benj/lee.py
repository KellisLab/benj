### Use squidpy generic neighbors
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix, diags, issparse, isspmatrix_csr
from scipy.sparse.linalg import eigsh
from spatialdata import SpatialData
from tqdm.auto import tqdm

from squidpy._constants._pkg_constants import Key
from squidpy._docs import d, inject_docs
from squidpy.gr._utils import _assert_connectivity_key, _assert_non_empty_sequence


def lee(
    adata: AnnData | SpatialData,
    n_pcs_pos: int = 25,
    n_pcs_neg: int = 25,
    genes: str | Sequence[str] | None = None,
    key_added: str = "X_lee",
    connectivity_key: str = Key.obsp.spatial_conn(),
    zero_center: bool = True,
    scale: bool = True,
    self_neighbor_weight: bool = False,
    batch_size: int = 1000,
) -> None:
    """Compute Lee's L spatially weighted PCA

    Lee's L is a weighted covariance structure based on spatially lagged covariance, similar to Moran's I, but can be used as a dimensionality reduction technique in place of PCA.
    See :cite:`lee2001` for reference.

    Parameters
    ----------
    %(adata)s
    n_pcs_pos
        Number of positive principal components to take. The covariance can have some strange structure so it can make sense to include negative PCs.
    n_neg_pcs
        Number of negative principal components to take. The covariance can have some strange structure so it can make sense to include negative PCs.
    genes
        A mask variable similar to squidpy.gr.spatial_autocorr. TODO implement
    key_added
        Where to store the embeddings and loadings
    connectivity_key
    
    Returns
    -------
    None
    """
    if isinstance(adata, SpatialData):
        adata = adata.table
    _assert_connectivity_key(adata, connectivity_key)
    if genes is None:
        genes = adata.var_names.values
    else:
        raise NotImplementedError()
    genes = _assert_non_empty_sequence(genes, name="genes")
    W = adata.obsp[connectivity_key]
    if not isspmatrix_csr(W):
        W = csr_matrix(W)
    if self_neighbor_weight:
        W.setdiag(np.ones(W.shape[0], dtype=W.dtype))
    W.eliminate_zeros()
    W = W.T.dot(W) * (W.shape[0] / np.ravel(W.sum(0) ** 2).sum())  ### Denominator is 1(V'V)1=squared norm
    ### extract means
    mu = np.ravel(adata.X.mean(0))
    cov = np.zeros((adata.shape[1], adata.shape[1]), dtype=np.float32)
    zero_center |= scale
    inv_sd = np.zeros(adata.shape[1]) * np.nan
    for ileft in tqdm(np.arange(0, len(all_indices), batch_size), desc="Computing covariance"):
        iright = min(ileft + batch_size, adata.shape[1])
        X = adata.X[:, ileft:iright].copy()
        if issparse(X):
            X = X.todense().A
        X = np.asarray(X)
        if zero_center:
            X = X - mu[None, ileft:iright]
        if scale:
            sd = np.std(X, axis=0).clip(1e-25, np.inf)
            inv_sd[ileft:iright] = 1 / sd
            X = X * inv_sd[None, ileft:iright]
        XW = W.dot(X).T
        for jleft in tqdm(np.arange(0, iright, batch_size), leave=False):
            jright = min(jleft + batch_size, iright)
            Y = adata.X[:, jleft:jright].copy()
            if issparse(Y):
                Y = Y.todense().A
            Y = np.asarray(Y)
            if zero_center:
                Y = Y - mu[None, jleft:jright]
            if scale:
                Y = Y * inv_sd[None, jleft:jright]
            cov[ileft:iright, jleft:jright] = XW.dot(Y)
            if jleft < ileft:  ### don't recompute if not the same block
                cov[jleft:jright, ileft:iright] = cov[ileft:iright, jleft:jright].T
    if n_pcs_pos == 0:
        e_val, e_vec = eigsh(cov / (adata.shape[0] * 2 - 2), n_pcs_neg, which="SA")
    elif n_pcs_neg == 0:
        e_val, e_vec = eigsh(cov / (adata.shape[0] * 2 - 2), n_pcs_pos, which="LA")
    else:
        e_val, e_vec = eigsh(cov / (adata.shape[0] * 2 - 2), 2 * max(n_pcs_pos, n_pcs_neg), which="BE")
        if n_pcs_pos != n_pcs_neg:
            idx = np.argsort(e_val)
            idx = np.hstack((idx[:n_pcs_neg], idx[-n_pcs_pos:]))
            e_val, e_vec = e_val[idx], e_vec[:, idx]
    adata.varm[key_added] = e_vec
    adata.uns["lee"] = {
        "eigenvalues": e_val,
        "zero_center": zero_center,
        "scale": scale,
        "n_pcs_pos": n_pcs_pos,
        "n_pcs_neg": n_pcs_neg,
        "key_added": key_added,
        "connectivities_key": connectivity_key,
    }
    if scale:
        ### Pre-divide by using eigenvectors as dim is smaller
        e_vec = diags(inv_sd).dot(e_vec)
    adata.obsm[key_added] = adata.X.dot(e_vec)
    if zero_center:
        adata.obsm[key_added] -= mu[None, :] @ e_vec
