
class WeightedCCA:
    def __init__(self, n_components:int=2, z_scale:bool=True, normalize_weights:bool=True):
        from sklearn.preprocessing import StandardScaler
        from sklearn.cross_decomposition import CCA
        self.z_scale = z_scale
        self.normalize_weights = normalize_weights
        self.cca = CCA(n_components=n_components)
        self.scaler_X, self.scaler_Y = StandardScaler(), StandardScaler()
    def _normalize_weights(self, W):
        import numpy as np
        W = np.asarray(W)
        W[np.isnan(W)] = 0
        if self.normalize_weights:
            return W / np.linalg.norm(W, 2)
        return W
    def fit(self, X, Y, W=None):
        import numpy as np
        if W is None:
            W = np.ones(X.shape[0])
        W = self._normalize_weights(W)
        if self.z_scale:
            X = self.scaler_X.fit_transform(X)
            Y = self.scaler_Y.fit_transform(Y)
        X_w = X * np.sqrt(W).reshape(-1, 1)
        Y_w = Y * np.sqrt(W).reshape(-1, 1)
        self.cca.fit(X_w, Y_w)
        return self
    def transform(self, X, Y, W=None):
        import numpy as np
        if W is None:
            W = np.ones(X.shape[0])
        W = self._normalize_weights(W)
        if self.z_scale:
            X = self.scaler_X.transform(X)
            Y = self.scaler_Y.transform(Y)
        X_w = X * np.sqrt(W).reshape(-1, 1)
        Y_w = Y * np.sqrt(W).reshape(-1, 1)
        return self.cca.transform(X_w, Y_w)
    def fit_transform(self, X, Y, W=None):
        return self.fit(X, Y, W).transform(X, Y, W)
