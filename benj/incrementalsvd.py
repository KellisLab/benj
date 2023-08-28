
class IncrementalSVD:
        """Strategy:
           1. Calculate residuals of new batch after projecting singular vectors.
           2. QR the residuals to get orthogonal vectors
           3. Create augmented matrix by placing the old singular values on the diagonal and appending orthogonalized new batch.
           4. Compute the SVD of the augmented matrix.
           5. Update singular values and vectors.
        """
        def __init__(self, n_components:int=50):
                self.n_components = n_components
                self.s = None
                self.V = None
        def partial_fit(self, X, method_qr:bool=False):
                import numpy as np
                import scipy.sparse
                from scipy.sparse.linalg import svds
                if self.s is None and self.V is None:
                        U, s, VT = svds(X, k=self.n_components)
                        _, self.s, self.V = U[:, ::-1], s[::-1], VT[::-1, :].T
                elif method_qr:
                        residuals = X - np.dot(np.dot(X, self.V), self.V.T)
                        Q, R = np.linalg.qr(residuals)
                        del residuals
                        aug_upper = np.hstack((np.diag(self.s),
                                               np.zeros((len(self.s), R.shape[1] - len(self.s)))))
                        aug_lower = R
                        augmented = np.vstack((aug_upper, aug_lower))
                        del aug_upper, aug_lower
                        Uz, sz, VzT = svds(augmented, k=self.n_components)
                        self.s = sz[::-1]
                        self.V = VzT[::-1, :].T
                else:
                        X = scipy.sparse.vstack((np.diag(self.s).dot(self.V.T),
                                                 X))
                        Uz, sz, VzT = svds(X, k=self.n_components)
                        self.s = sz[::-1]
                        self.V = VzT[::-1, :].T
        def transform(self, X):
                import numpy as np
                return np.dot(X, self.V)
