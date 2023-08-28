
from typing import Union
class IncrementalLSI:
        def __init__(self, var_sums,
                     n_comps:int=50,
                     log_tf:bool=True, log_idf:bool=True, log_tfidf:bool=False,
                     scale_factor:Union[int, float]=1e4,
                     ):
                import numpy as np
                from .incrementalsvd import IncrementalSVD
                if log_tfidf and (log_tf or log_idf):
                        raise AttributeError(
                                "When returning log(TF*IDF), \
                                applying neither log(TF) nor log(IDF) is possible."
                        )
                var_sums = np.ravel(var_sums)
                self.tf = np.zeros(len(var_sums), dtype="f8")
                self.tf = np.divide(scale_factor, var_sums, where=var_sums > 0, out=self.tf)
                self.log_tf = log_tf
                self.log_idf = log_idf
                self.log_tfidf = log_tfidf
                self.svd = IncrementalSVD(n_comps)
        def partial_fit(self, X):
                self.svd.partial_fit(self._tfidf(X))
                return self
        def transform(self, X):
                return self.svd.transform(self._tfidf(X))
        def _tfidf(self, X):
                import numpy as np
                from scipy.sparse import issparse, dia_matrix, csr_matrix
                idf = np.ravel(X.shape[0] / X.sum(0))
                if self.log_idf:
                        idf = np.log1p(idf)
                if issparse(X):
                        tf = np.dot(dia_matrix((self.tf, 0),
                                               shape=(self.tf.shape[0],
                                                      self.tf.shape[0])),
                                    X)
                else:
                        tf = X / self.tf
                if self.log_tf:
                        tf = np.log1p(tf)
                if issparse(tf):
                        idf = dia_matrix((idf, 0), shape=(idf.size, idf.size))
                        tf_idf = np.dot(tf, idf)
                else:
                        tf_idf = np.dot(csr_matrix(tf), csr_matrix(np.diag(idf)))
                if self.log_tfidf:
                        tf_idf = np.log1p(tf_idf)
                return np.nan_to_num(tf_idf, 0)
