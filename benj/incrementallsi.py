
from typing import Union
class IncrementalTFIDF:
        def __init__(self, var_means,
                     log_tf:bool=True, log_idf:bool=True, log_tfidf:bool=False,
                     scale_factor:Union[int, float]=1e4):
                import numpy as np
                if log_tfidf and (log_tf or log_idf):
                        raise AttributeError(
                                "When returning log(TF*IDF), \
                                applying neither log(TF) nor log(IDF) is possible."
                        )
                var_means = np.ravel(var_means)
                self.idf = np.zeros(len(var_means), dtype="f8")
                self.idf = np.divide(1, var_means, where=var_means > 0, out=self.idf)
                self.scale_factor = scale_factor
                self.log_tf = log_tf
                if log_idf:
                        self.idf = np.log1p(self.idf)
                self.log_tfidf = log_tfidf
        def transform(self, X):
                import numpy as np
                from scipy.sparse import issparse, dia_matrix, csr_matrix
                ### First get X sums per row
                xs1 = np.ravel(X.sum(1))
                tf = np.zeros(len(xs1), dtype="f8")
                tf = np.divide(self.scale_factor, xs1, where=xs1>0, out=tf)
                del xs1
                if issparse(X):
                        tf = np.dot(dia_matrix((tf, 0),
                                               shape=(tf.shape[0],
                                                      tf.shape[0])),
                                    X)
                else:
                        tf = np.reshape(tf, (-1, 1))
                        tf = X * tf
                if self.log_tf:
                        tf = np.log1p(tf)
                if issparse(tf):
                        idf = dia_matrix((self.idf, 0), shape=(self.idf.size, self.idf.size))
                        tf_idf = np.dot(tf, idf)
                else:
                        tf_idf = np.dot(csr_matrix(tf), csr_matrix(np.diag(self.idf)))
                if self.log_tfidf:
                        tf_idf = np.log1p(tf_idf)
                return np.nan_to_num(tf_idf, 0)

class IncrementalLSI:
        def __init__(self, var_means,
                     n_comps:int=50,
                     log_tf:bool=True, log_idf:bool=True, log_tfidf:bool=False,
                     scale_factor:Union[int, float]=1e4,
                     ):
                from .incrementalsvd import IncrementalSVD
                self.tfidf = IncrementalTFIDF(var_means=var_means,
                                              log_tf=log_tf, log_idf=log_idf,
                                              log_tfidf=log_tfidf,
                                              scale_factor=scale_factor)
                self.svd = IncrementalSVD(n_comps)
        def partial_fit(self, X):
                self.svd.partial_fit(self.tfidf.transform(X))
                return self
        def transform(self, X):
                return self.svd.transform(self.tfidf.transform(X))
