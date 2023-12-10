
from typing import Union, List, Optional
from pathlib import Path
_PathLike=Union[str, Path]

def aggregate_collection(adata, which:Union[str, List[str]]="X", view:bool=True):
    from tqdm.auto import tqdm
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import anndata
    from anndata.experimental.multi_files import AnnCollection
    if isinstance(which, str):
        which = [which]
    tbl = {}
    good_samples = pd.unique(adata.obs[adata.uns["H5AD"]["sample_key"]])
    for k, v in tqdm(adata.uns["H5AD"]["files"].items()):
        if k not in good_samples:
                continue
        tbl[k] = anndata.read_h5ad(v, backed="r")
        del tbl[k].raw
        if "layers" not in which and "all" not in which:
                del tbl[k].layers
    ac = AnnCollection(tbl, join_vars="inner")
    if view:
        return ac[adata.obs_names, adata.var_names]
    else:
        ac


def aggregate_var(tbl:dict):
    import re
    import numpy as np
    import pandas as pd
    import anndata
    def _aggregate_stats(df_tbl, prefix=""):
        ### get number of cells per .var
        nf = {s: df.get("total_ncells", df.get("n_cells_by_counts")).max() for s, df in df_tbl.items()}
        mf = pd.concat({s: df[prefix + "mean"] for s, df in df_tbl.items()}, axis=1)
        sf = pd.concat({s: df[prefix + "std"] for s, df in df_tbl.items()}, axis=1)
        sf = sf.loc[mf.index.values, mf.columns.values]
        nf = pd.Series(nf, index=mf.columns.values)
        overall_mean = (mf * nf).sum(1) / nf.sum()
        weighted_var = (nf - 1) * sf * sf
        mean_diff_sq = nf * (mf - overall_mean.values[:, None])**2
        overall_var = (weighted_var + mean_diff_sq).sum(1) / (nf.sum() - len(df_tbl))
        return pd.DataFrame({prefix + "mean": overall_mean, prefix + "std": np.sqrt(overall_var)})
    def _aggregate_hvg(gb,
                       min_disp: Optional[float] = 0.5,
                       max_disp: Optional[float] = np.inf,
                       min_mean: Optional[float] = 0.0125,
                       max_mean: Optional[float] = 3,):
        hf = gb.agg({
            "means": np.nanmean,
            "dispersions": np.nanmean,
            "dispersions_norm": np.nanmean,
            "highly_variable": np.nansum})
        hf.rename(columns={"highly_variable": "highly_variable_nbatches"}, inplace=True)
        dispersion_norm = hf.dispersions_norm.values
        dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
        gene_subset = np.logical_and.reduce(
            (
                hf.means > min_mean,
                hf.means < max_mean,
                hf.dispersions_norm > min_disp,
                hf.dispersions_norm < max_disp,
            )
        )
        hf['highly_variable'] = gene_subset
        return hf
    comm_cols = None
    var_tbl = {}
    for k, data in tbl.items():
        if isinstance(data, anndata.AnnData):
            var_tbl[k] = data.var.copy()
            var_tbl[k]["total_ncells"] = data.shape[0]
        elif isinstance(data, pd.DataFrame):
            var_tbl[k] = data.copy()
        else:
            continue
        var_tbl[k]["gene"] = var_tbl[k].index.values
    cf = pd.concat(var_tbl, axis=0)
    var = cf.groupby("gene").agg({"total_counts": np.nansum,
                                  "n_cells_by_counts": np.nansum})
    var["log1p_total_counts"] = np.log1p(var["total_counts"])
    if "dispersions" in cf.columns:
        dvar = _aggregate_hvg(cf.groupby("gene"))
        for col in dvar.columns:
            var[col] = dvar.loc[var.index.values, col].values
    if "std" in cf.columns and "mean" in cf.columns:
        svar = _aggregate_stats(var_tbl)
        for col in svar.columns:
            var[col] = svar.loc[var.index.values, col].values
    layers = set([re.sub("_mean$", "", s) for s in cf.columns if s.endswith("_mean")])
    layers &= set([re.sub("_std$", "", s) for s in cf.columns if s.endswith("_std")])
    for lay in layers:
        try:
            svar = _aggregate_stats(var_tbl, prefix="%s_" % lay)
            for col in svar.columns:
                var[col] = svar.loc[var.index.values, col].values
        except:
            print("layer:", lay)
            print("mean:", [k for k, var in var_tbl.items() if "%s_mean" % lay not in var.columns ])
            print("std:", [k for k, var in var_tbl.items() if "%s_std" % lay not in var.columns ])
            pass
    return var

def aggregate_load(adata, which:Union[str, List[str]]="X"):
    import scanpy as sc
    from .timer import template as stopwatch
    sw = stopwatch()
    if isinstance(which, str):
        which = [which]
    ac = aggregate_collection(adata, which=which)
    if "X" in which or "all" in which:
        with sw("Loading X (%d, %d)" % adata.shape):
            adata.X = ac.X
    if "layers" in which or "all" in which:
        with sw("Loading layers (%d, %d)" % adata.shape):
            adata.layers = ac.layers
    return adata

def aggregate_concat(metadata=None, directory:Union[_PathLike, List[_PathLike]]=None,
                     h5ad:Union[_PathLike, List[_PathLike]]=None,
                     sample_key="Sample", calc_qc:bool=True,
                     min_cells_per_sample:int=30,
                     sep="\t", verbose:bool=True,
                     **kwargs):
    """Metadata+directory, or h5ad with or without metadata"""
    import os
    import numpy as np
    import pandas as pd
    from tqdm.auto import tqdm
    import anndata
    from .timer import template as stopwatch
    from .load_anndata import find_sample
    sw = stopwatch()
    if h5ad is not None:
        if isinstance(h5ad, str) or isinstance(h5ad, Path):
            h5ad = [h5ad]
        h5ad = {os.path.basename(p).replace(".h5ad", ""): p for p in h5ad}
    else:
        h5ad = {}
    if metadata is None:
        ### options: pass H5AD with or without metadata,
        ### or pass directory with metadata
        metadata = pd.DataFrame(index=list(h5ad.keys()))
    elif not isinstance(metadata, pd.DataFrame):
        if os.path.exists(metadata):
            metadata = pd.read_csv(metadata, sep=sep, index_col=0, low_memory=False)
        else:
            raise ValueError("Metadata does not exist!")
    if directory is None and len(h5ad) == 0:
        raise ValueError("Must provide either directory or h5ad")
    if verbose and isinstance(metadata, pd.DataFrame):
        print("Metadata:")
        print(metadata)
    adata_tbl = {}
    fname_tbl = {}
    scrub_tbl = {}
    h5ad_tbl = {}
    bad = []
    with sw("Reading H5AD files"):
        for sample in tqdm(metadata.index.values):
            ### TODO sample may be aggr not actual sample
            ### find sample checks 1: h5ad 2: directory+metadata, 3: sample.h5ad
            adata, fname = find_sample(directory=directory, h5ad=h5ad.get(sample),
                                       sample=sample, metadata=metadata, qc=calc_qc,
                                       min_cells_per_sample=min_cells_per_sample,
                                       verbose=verbose,
                                       **kwargs)
            if adata is None:
                bad.append(sample)
            else:
                adata_tbl[sample] = adata
                fname_tbl[sample] = fname
                if "H5AD" in adata.uns and adata.uns["H5AD"].get("sample_key", "") == sample_key:
                    ### Aggregating raw
                    h5ad_tbl |= adata.uns["H5AD"]["files"]
                else: ### already set, using wrong sample/md combo
                    adata.obs[sample_key] = sample
                if "scrublet" in adata.uns:
                    if "batches" in adata.uns["scrublet"]:
                        scrub_tbl |= adata.uns["scrublet"]["batches"]
                    else:
                        scrub_tbl[sample] = adata.uns["scrublet"]
    if bad:
        print("Bad samples: ", ",".join(bad))
    total_cells = np.sum([adata.shape[0] for _, adata in adata_tbl.items()])
    with sw("Concatenating %d cells into one AnnData object" % total_cells):
        if calc_qc:
            calc_qc = aggregate_var(adata_tbl)
        adata = anndata.concat(adata_tbl, merge="same", uns_merge="same")
        tk = list(adata_tbl.keys())
        del adata_tbl
        if len(tk) == len(scrub_tbl.keys()):
            adata.uns["scrublet"] = {"batches": scrub_tbl,
                                     "batched_by": sample_key}
        if len(h5ad_tbl) > 0:
            adata.uns["H5AD"] = {"sample_key": sample_key, "files": h5ad_tbl}
        if isinstance(calc_qc, pd.DataFrame):
            for col in calc_qc.columns:
                adata.var[col] = calc_qc.loc[adata.var_names, col].values
    if "H5AD" not in adata.uns:
        ### Creating raw
        adata.uns["H5AD"] = {"sample_key": sample_key,
                             "files": fname_tbl}
    return adata
