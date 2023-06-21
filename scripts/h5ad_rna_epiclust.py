#!/usr/bin/env python3

def run(adata, output:str=None, batch:str=None, hvg:int=5000, resolution:float=1., pc:bool=True,
        max_comm_size:int=2500, min_comm_size:int=3, min_cells:int=3, n_neighbors:int=10,
        umap:bool=True, filter_z:float=2., **kwargs):
        import numpy as np
        import scanpy as sc
        import epiclust as ec
        if pc and "pc" in adata.var.columns:
                adata = adata[:, adata.var["pc"]]
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=hvg, subset=True, batch_key=batch)
        if batch is not None and batch in adata.obs.columns:
                sc.pp.combat(adata)
        else:
                sc.pp.scale(adata)
        sc.pp.pca(adata, zero_center=False)
        for power in np.linspace(0, 1, 5):
                print("power:", power)
                ec.fit(adata, power=power)
                ec.neighbors(adata, n_neighbors=n_neighbors)
        ec.filter_var(adata, adata.uns["epiclust"]["graphs"], z=filter_z, min_cells=min_cells)
        ec_kwargs = {}
        if max_comm_size > 0:
                ec_kwargs["max_comm_size"] = max_comm_size
        ec.leiden(adata, adata.uns["epiclust"]["graphs"],
                  min_comm_size=max(min_comm_size, 0),
                  resolution=resolution, **ec_kwargs)
        if umap:
                ec.umap(adata, adata.uns["epiclust"]["graphs"])
                adata.var["UMAP0"] = adata.varm["X_umap"][:, 0]
                adata.var["UMAP1"] = adata.varm["X_umap"][:, 1]
        adata.var.loc[~adata.var["leiden"].isna(), :].to_csv(args["output"], sep="\t")

if __name__ == "__main__":
    import argparse
    import benj
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-i", "--input", dest="h5ad", required=True)
    ap.add_argument("--hvg", type=int, default=5000)
    ap.add_argument("-n", "--neighbors", type=int, default=10)
    ap.add_argument("-b", "--batch")
    ap.add_argument("-r", "--resolution", type=float, default=1.)
    ap.add_argument("--max-comm-size", type=int, default=2500)
    ap.add_argument("--min-comm-size", type=int, default=3)
    ap.add_argument("--min-cells", type=int, default=3)
    ap.add_argument("--pc", dest="pc", action="store_true")
    ap.add_argument("--all-genes", dest="pc", action="store_false")
    ap.add_argument("--umap", dest="umap", action="store_true")
    ap.add_argument("--no-umap", dest="umap", action="store_false")
    ap.add_argument("-z", "--filter-z", type=float, default=2.)
    ap.set_defaults(pc=True, umap=True)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    sw = benj.stopwatch()
    with sw("Reading H5AD"):
        adata = benj.parse_anndata(**args)
    run(adata, **args)
