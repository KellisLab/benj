#!/usr/bin/env python3

def h5ad_per_sample_to_xlsx(adata, output, target_sum:float=None, **kwargs):
    """
    sheet 1: obs (pseudobulk)
    sheet 2: counts (pseudobulk)
    sheet 3: quantile 0
    sheet 4: quantile 0.25
    sheet 5: quantile 0.5
    sheet 6: quantile 0.75
    sheet 7: quantile 1
    """
    import numpy as np
    from tqdm.auto import tqdm
    import pandas as pd
    import scipy.sparse
    import scanpy as sc
    from benj import pseudobulk
    import benj
    sw = benj.stopwatch()
    with sw("Pseudobulking"):
            pdata = pseudobulk(adata, ["Sample"])
            pdata.obs.index = pdata.obs["Sample"].values
            qc_vars = set(["mt", "ribo", "pc", "chrX", "chrY"]) & set(pdata.var.columns)
            sc.pp.calculate_qc_metrics(pdata, qc_vars=qc_vars, inplace=True)
            pdata.layers["raw"] = pdata.X.copy()
            sc.pp.normalize_total(pdata)
            sc.pp.log1p(pdata)
            sc.pp.scale(pdata)
            sc.pp.pca(pdata, zero_center=False)
            pdata.X = pdata.layers["raw"].copy()
            del pdata.layers["raw"]
            pdata.var.index.name = "gene"
    q = pd.Series([0, 0.25, 0.5, 0.75, 1])
    q.index = ["Quantile %.2f" % x for x in q.values]
    P = np.zeros((pdata.shape[0], len(q), pdata.shape[1]))
    with sw("Calculating quantiles"):
        for i in tqdm(np.arange(pdata.shape[0])):
            sample = pdata.obs.index.values[i]
            xdata = adata[adata.obs["Sample"] == sample, :].copy()
            if scipy.sparse.issparse(xdata.X):
                xdata.X = np.asarray(xdata.X.todense())
            sc.pp.normalize_total(xdata, target_sum=target_sum)
            sc.pp.log1p(xdata)
            P[i, :, :] = np.quantile(xdata.X, q, axis=0)
            IQR = P[i, 3, :] - P[i, 1, :]
            P[i, 0, :] = P[i, 1, :] - 1.5 * IQR
            P[i, 4, :] = P[i, 3, :] + 1.5 * IQR
            ### Get masked max, min to remove 'outliers' from boxplot
            ma = np.ma.masked_array(xdata.X, xdata.X > P[i, [0], :]).min(0)
            P[i, 0, ~ma.mask] = ma.compressed()
            ma = np.ma.masked_array(xdata.X, xdata.X < P[i, [4], :]).max(0)
            P[i, 4, ~ma.mask] = ma.compressed()
            del xdata
    with sw("Writing output"):
        with pd.ExcelWriter(output) as writer:
            pf = pd.DataFrame(pdata.obsm["X_pca"],
                                          index=pdata.obs_names,
                                          columns=["PC%d" % i for i in range(pdata.obsm["X_pca"].shape[1])])
            obs = pd.concat((pdata.obs, pf), axis=1)
            obs.to_excel(writer, sheet_name="sample_metadata")
            pdata.var.to_excel(writer, sheet_name="gene_info")
            pd.DataFrame(pdata.X.T,
                         columns=pdata.obs_names,
                         index=pdata.var_names).to_excel(writer, sheet_name="counts")
            for i, quant in enumerate(q):
                pd.DataFrame(P[:, i, :].T,
                             columns=pdata.obs_names,
                             index=pdata.var_names).to_excel(writer,
                                                             sheet_name=q.index.values[i])

if __name__ == "__main__":
    import argparse
    import benj
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", dest="h5ad", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    adata = benj.parse_anndata(**args)
    h5ad_per_sample_to_xlsx(adata, **args)
