
from typing import Union, List
from pathlib import Path
_PathLike=Union[str, Path]

def load_anndata(h5ad:_PathLike,
                 obs_annotation:Union[_PathLike, List[_PathLike]]=None,
                 subset:dict=None, qc:bool=True,
                 elt:Union[str, List[str]]=None, logger="benj", verbose:bool=True,
                 obs_min:dict=None, obs_max:dict=None, sep:str="\t", **kwargs):
    import numpy as np
    import pandas as pd
    import anndata
    from .utils import read_elems
    if isinstance(obs_annotation, str):
        obs_annotation = pd.read_csv(obs_annotation, sep=sep, index_col=0, low_memory=False)
    elif isinstance(obs_annotation, list):
        obs_annotation = pd.concat([pd.read_csv(fname, sep=sep, index_col=0, low_memory=False) for fname in obs_annotation], axis=1, join="inner")
    if elt is None:
        elt = ["obs", "var", "obsm", "varm", "varp", "uns"]
    elif elt == "all":
        elt = ["obs", "var", "obsm", "varm", "varp", "uns", "X", "layers", "raw"]
    temp_X = qc and "X" not in elt
    if temp_X:
        elt.append("X")
    adata = anndata.AnnData(**read_elems(h5ad, elt))
    obs = adata.obs.copy()
    if isinstance(obs_annotation, pd.DataFrame):
        obs_annotation = obs_annotation.loc[obs_annotation.index.isin(obs.index.values), :]
        obs = obs.loc[obs_annotation.index.values, :]
        for cn in obs_annotation.columns.values:
            obs[cn] = obs_annotation[cn].values
    if isinstance(subset, dict):
        for sub_k, sub_v in subset.items():
            if not isinstance(sub_v, list):
                sub_v = [sub_v]
            dt = obs[sub_k].dtype
            if dt.name == "category":
                dt = str
            sub_v = np.asarray(sub_v, dt)
            obs = obs.loc[obs[sub_k].isin(sub_v), :]
            if obs.shape[0] == 0:
                raise RuntimeError("%s not in %s" % (",".join(sub_v), sub_k))
    if isinstance(obs_min, dict):
        for min_k, min_v in obs_min.items():
            if min_k in obs.columns:
                obs = obs.loc[obs[min_k] >= min_v, :]
    if isinstance(obs_max, dict):
        for max_k, max_v in obs_max.items():
            if max_k in obs.columns:
                obs = obs.loc[obs[max_k] <= max_v, :]
    adata = adata[obs.index.values, :].copy()
    for cn in set(obs.columns) - set(adata.obs.columns):
        adata.obs[cn] = obs[cn]
    if qc:
        import scanpy as sc
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        if "gene_ids" in adata.var.columns and adata.shape[0] > 0:
            import anndata
            bdata = anndata.AnnData(sc.pp.normalize_total(adata, target_sum=10000, inplace=False)["X"],
                                    obs=adata.obs, var=adata.var)
            sc.pp.log1p(bdata)
            sc.pp.highly_variable_genes(bdata)
            sc.pp.scale(bdata)
            adata.var = bdata.var
            del bdata
    if temp_X:
        del adata.X
    if verbose:
        import logging
        logger = logging.getLogger(logger)
        logging.basicConfig(level=logging.INFO)
        logger.info("Read %d cells" % adata.shape[0])
    return adata

def find_sample(sample:str, metadata=None, directory:Union[_PathLike, List[_PathLike]]=None,
                h5ad:_PathLike=None, qc:bool=True,
                min_cells_per_sample:int=30,
                **kwargs):
    import os
    import anndata
    from .load_anndata import load_anndata
    if h5ad is not None:
        fname = h5ad
    elif directory is not None:
        if not isinstance(directory, list):
            directory = [directory]
        for dname in directory:
            fname = os.path.join(dname, "%s.h5ad" % sample)
            if os.path.isfile(fname):
                break
        if not os.path.isfile(fname):
            raise RuntimeError("Sample %s.h5ad did not exist in the directories: %s" % (sample, ",".join(directory)))
    else:
        fname = "%s.h5ad" % sample
    adata = load_anndata(h5ad=fname, **kwargs)
    if metadata is not None:
        for cn in metadata.columns:
            adata.obs[cn] = metadata.loc[sample, cn]
    if adata.shape[0] >= min_cells_per_sample:
        return adata, os.path.abspath(fname)
    else:
        return None, fname
