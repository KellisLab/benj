class _RawNormalizeLogScaler:
        def __init__(self, mean, std, target_sum:int=None):
                self.mean = mean
                self.std = std
                self.target_sum = target_sum
        def transform(self, data):
                import numpy as np
                import pandas as pd
                import anndata
                import scanpy as sc
                if not isinstance(data, anndata.AnnData):
                        data = anndata.AnnData(data)
                X = sc.pp.normalize_total(data, target_sum=self.target_sum, inplace=False)["X"]
                X = sc.pp.log1p(X)
                X = (X - np.reshape(self.mean, (1, -1))) / np.reshape(self.std, (1, -1))
                return np.asarray(X)

class IncrementalPCA:
        def __init__(self, mean, std, target_sum:int=None, n_components:int=50):
                from .incrementalsvd import IncrementalSVD
                self.scaler = _RawNormalizeLogScaler(mean=mean, std=std, target_sum=target_sum)
                self.svd = IncrementalSVD(n_components=n_components)
        def partial_fit(self, X):
                self.svd.partial_fit(self.scaler.transform(X))
                return self
        def transform(self, X):
                return self.svd.transform(self.scaler.transform(X))
