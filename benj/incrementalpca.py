class _RawNormalizeLogScaler:
        def __init__(self, mean, std, target_sum:int=None, eps=1e-50, pre_normalized:bool=False):
                import numpy as np
                self.mean = np.ravel(mean).reshape(1, -1)
                self.inv_std = np.ravel(1. / std).reshape(1, -1).clip(eps, np.inf)
                self.target_sum = target_sum
                self.pre_normalized = pre_normalized
                if target_sum is None: ### weird python thing
                        self.target_sum = None
        def transform(self, data):
                import numpy as np
                import pandas as pd
                import anndata
                import scanpy as sc
                from scipy.sparse import issparse
                if not isinstance(data, anndata.AnnData):
                        data = anndata.AnnData(data)
                if self.pre_normalized:
                        X = data.X
                else:
                        X = sc.pp.normalize_total(data, target_sum=self.target_sum, inplace=False)["X"]
                        X = sc.pp.log1p(X)
                if issparse(X):
                        X = (X.toarray() - self.mean) * self.inv_std
                else:
                        X = (X - self.mean) * self.inv_std 
                return np.asarray(X)

class IncrementalPCA:
        def __init__(self, var, target_sum:int=None, n_components:int=50, pre_normalized:bool=False):
                from .incrementalsvd import IncrementalSVD
                self.scaler = _RawNormalizeLogScaler(mean=var["mean"], std=var["std"], target_sum=target_sum, pre_normalized=pre_normalized)
                self.svd = IncrementalSVD(n_components=n_components)
        def partial_fit(self, X):
                import gc
                self.svd.partial_fit(self.scaler.transform(X))
                gc.collect()
                return self
        def transform(self, X):
                import gc
                out = self.svd.transform(self.scaler.transform(X))
                gc.collect()
                return out
