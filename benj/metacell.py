
def resample(obs, batch_key:str=None, min_batch_count:int=10, seed:int=None):
    import numpy as np
    if seed is not None:
            np.random.seed(seed)
    I = np.unique(np.random.choice(obs.index, obs.shape[0], replace=True))
    if batch_key in obs.columns:
        obs = obs.loc[I, [batch_key]]
        vc = obs[batch_key].value_counts()
        vc = vc[vc >= min_batch_count]
        obs = obs.loc[obs[batch_key].isin(vc.index),:]
    return obs.index

def over_cluster_resolution(resolution, n_obs):
    if resolution is None:
        if n_obs < 5000:
            resolution = 5
        elif n_obs < 20000:
            resolution = 10
        elif n_obs < 40000:
            resolution = 15
        elif n_obs < 100000:
            resolution = 20
        elif n_obs < 200000:
            resolution = 25
        else:
            resolution = 30
    return resolution

def _mc_prepare(adata, idx, neighbors_key:str=None):
        import scanpy as sc
        import anndata
        flag = adata.obs_names.isin(idx)
        neighbors_key = "neighbors" if neighbors_key is None else neighbors_key
        P = adata.uns[neighbors_key].get("params")
        rep = P.get("use_rep", "X_pca")
        xdata = anndata.AnnData(obsm={rep: adata.obsm[rep][flag,:]}, obs=adata.obs.iloc[flag,:])
        method = P.get("method", "umap")
        metric = P.get("metric", "euclidean")
        n_neighbors = P.get("n_neighbors", 15)
        n_pcs = P.get("n_pcs", 50)
        sc.pp.neighbors(xdata, n_neighbors=n_neighbors,
                        n_pcs=n_pcs, use_rep=rep,
                        metric=metric, method=method)
        return xdata

def metacell_mudata(mdata, idx, resolution:float=None, neighbors_key:str=None,
                     max_comm_size:int=0, n_leiden_iterations:int=-1, layer_weights=None):
        from .utils import leiden_multiplex
        import muon as mu
        xdata = {}
        for k, adata in mdata.mod.items():
                xdata[k] = _mc_prepare(adata, idx=idx, neighbors_key=neighbors_key)
        xdata = mu.MuData(xdata)
        if layer_weights is not None:
                if isinstance(layer_weights, dict):
                        layer_weights = [layer_weights[k] for k in xdata.mod.keys()]
        leiden_multiplex(xdata, resolution=resolution,
                         max_comm_size=max_comm_size, n_iterations=n_leiden_iterations,
                         layer_weights=layer_weights)
        return xdata.obs["leiden"]

def metacell_anndata(adata, idx, resolution:float=None, neighbors_key:str=None,
                     max_comm_size:int=0, n_leiden_iterations:int=-1):
        import scanpy as sc
        xdata = _mc_prepare(adata, idx=idx, neighbors_key=neighbors_key)
        resolution = over_cluster_resolution(resolution=resolution, n_obs=xdata.shape[0])
        sc.tl.leiden(xdata, resolution=resolution, max_comm_size=max_comm_size, n_iterations=n_leiden_iterations)
        return xdata.obs["leiden"]

def metacell(data, niter:int=100,
             resolution:float=None,
             ### neighbors
             neighbors_key:str=None, ### requires neighbors() to be run beforehand
             ### leiden
             n_leiden_iterations:int=-1,
             max_comm_size:int=0,
             layer_weights=None):
        import pandas as pd
        import anndata
        from tqdm.auto import tqdm
         R = [resample(data.obs) for _ in range(niter)]
        out = []
        for idx in tqdm(R, desc="Computing metacells"):
                if isinstance(data, anndata.AnnData):
                        mc = metacell_anndata(data, idx, resolution=resolution, neighbors_key=neighbors_key, max_comm_size=max_comm_size, n_leiden_iterations=n_leiden_iterations)

                else:
                        import muon as mu
                        if isinstance(data, mu.MuData):
                                mc = metacell_mudata(data, idx, resolution=resolution, neighbors_key=neighbors_key, max_comm_size=max_comm_size, layer_weights=layer_weights, n_leiden_iterations=n_leiden_iterations)

                out.append(pd.get_dummies(mc).reindex(data.obs_names, fill_value=0))
        out = pd.concat(out, axis=1)
        out.columns = list(range(out.shape[1]))
        return out
