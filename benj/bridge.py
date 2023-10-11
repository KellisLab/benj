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
