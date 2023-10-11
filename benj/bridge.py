class Bridge:
        """usage: for multiome data, integrate each assay separately with single-omic
        then batch-correct.
        Pass dim_rep=batch corrected,
        and lap_rep as sklearn.manifold.spectral_embedding(mdata.obsp["connectivities"]) to lap_rep
        """
        def fit(self, dim_rep, lap_rep):
                import numpy as np
                self.w = np.linalg.pinv(dim_rep).dot(lap_rep)
                return self
        def transform(self, dim_rep):
                """Use unimodal (or all) """
                return dim_rep.dot(self.w)

def Laplacian(matrix, n=50, **kwargs):
        import numpy as np
        from scipy.sparse import diags
        from scipy.sparse.linalg import eigsh
        D_half_inv = diags(np.ravel(matrix.sum(1))**-0.5)
        L = -D_half_inv @ matrix @ D_half_inv
        L.setdiag(1 + L.diagonal())
        e_vals, e_vecs = eigsh(L, k=n + 1, which="SM", **kwargs)
        e_vals, e_vecs = e_vals[1:], e_vecs[:, 1:]
        return e_vals, e_vecs[:, ::-1]
