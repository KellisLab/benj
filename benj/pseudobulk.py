
def pseudobulk_valid_columns(obs, onehot):
        import pandas as pd
        if isinstance(onehot, pd.DataFrame):
                onehot = onehot.reindex(obs.index, fill_value=0)
        else:
                onehot = pd.get_dummies(obs.loc[:, onehot])
        goodcols = set(obs.columns.values)
        for i in tqdm(np.arange(onehot.shape[1]), desc="Checking pseudobulk cols"):
                idx_flag = onehot.iloc[:, i] > 0
                df = obs.loc[idx_flag, list(goodcols)]
                col_flag = df.apply(pd.Series.nunique, axis=0) == 1
                goodcols &= set(df.columns.values[col_flag])
        return list(goodcols)

def _pb_mul(batch, idx, mat, transform):
        import scipy.sparse
        import numpy as np
        import pandas as pd
        X = scipy.sparse.coo_matrix(batch)
        idx = np.asarray(idx)
        X = scipy.sparse.csr_matrix((X.data, (idx[X.row], X.col)), (transform.shape[0], batch.shape[1]))
        return transform.T.dot(X)

def _pb_anncollection(ac, transform, layer_list=None, batch_size:int=10000):
        import numpy as np
        import scipy.sparse
        import pandas as pd
        import anndata
        from tqdm.auto import tqdm
        from .utils import convert_X
        X = scipy.sparse.lil_matrix((transform.shape[1], ac.shape[1]), dtype="i8")
        layer_list = [] if layer_list is None else list(layer_list)
        L = {k: X.copy() for k in layer_list}
        for batch, idx in tqdm(ac.iterate_axis(batch_size),
                               total=np.ceil(ac.shape[0]/batch_size),
                               desc="Pseudobulking"):
                batch = batch[:, ac.var_names]
                X += _pb_mul(batch.X, idx=idx, mat=X, transform=transform)
                for lay in L.keys():
                        L[lay] += _pb_mul(batch.layers[lay], idx=idx, mat=L[lay], transform=transform)
        for lay in L.keys():
                L[lay] = convert_X(L[lay].tocsr())
        return anndata.AnnData(X=convert_X(X.tocsr()), layers=L, var=pd.DataFrame(index=ac.var_names))

def pseudobulk(adata, groupby, which:str="sum", batch_size:int=10000, prefix="PB#", dense:bool=False):
        import numpy as np
        import scipy.sparse
        import pandas as pd
        import anndata
        from anndata.experimental.multi_files import AnnCollection
        from .utils import convert_X
        if isinstance(groupby, pd.DataFrame):
                mc = groupby.reindex(adata.obs_names, fill_value=0)
        else:
                mc = pd.get_dummies(adata.obs.loc[:, groupby]).reindex(adata.obs_names, fill_value=0)
        cols = pseudobulk_valid_columns(adata.obs, mc)
        newobs = adata.obs.loc[mc.idxmax().values, cols]
        newobs.index = ["%s%d" % (prefix, x) for x, _ in enumerate(mc.columns)]
        if which == "sum":
                S = scipy.sparse.csr_matrix(mc.values, dtype="i8")
        elif which == "mean":
                S = scipy.sparse.csr_matrix(mc.values)
                S = scipy.sparse.diags(1/(mc.values.sum(1)+1e-100)).dot(S)
        del mc
        if isinstance(adata, AnnCollection):
                xdata = _pb_anncollection(adata, transform=S, layer_list=list(adata.layers.keys()), batch_size=batch_size)
        elif isinstance(adata, anndata.AnnData):
                if adata.isbacked:
                        ac = AnnCollection(adata)[adata.obs_names, adata.var_names]
                        xdata = _pb_anncollection(ac, transform=S, layer_list=list(adata.layers.keys()), batch_size=batch_size)
                else:
                        layers = {}
                        for lay in adata.layers.keys():
                                layers[lay] = convert_X(S.T.dot(adata.layers[lay]))
                        xdata = anndata.AnnData(X=convert_X(S.T.dot(adata.X)),
                                                var=adata.var,
                                                layers=layers)
        if dense:
                xdata.X = np.asarray(xdata.X.todense())
                for lay in xdata.layers.keys():
                        xdata.layers[lay] = np.asarray(xdata.layers[lay].todense())
        xdata.obs = newobs
        return xdata
