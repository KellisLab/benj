#!/usr/bin/env python3

def get_tss(tss:str):
    import pandas as pd
    tss = pd.read_csv(tsv, sep="\t", header=None)
    tss.columns = ["Chromosome", "Start", "End", "gene_id", "score", "strand"]
    df = tss.groupby(["Chromosome", "gene_id", "strand"]).agg(left=("Start", "min"),
                                                              right=("End", "max")).reset_index()
    df["interval"] = df["Chromosome"] + ":" + df["left"].astype(str) + "-" + df["right"].astype(str)
    df.index = df["gene_id"].values
    return df

def run(h5, output, sample:str=None, compression:int=9, tss:str=None, **kwargs):
    import os
    from warnings import warn
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import benj
    sw = benj.stopwatch()
    outs_dir = os.path.dirname(h5)
    bcfile = os.path.join(outs_dir, "filtered_feature_bc_matrix/barcodes.tsv.gz")
    vdir = os.path.join(outs_dir, "velocyto")
    vdata = None
    if not os.path.isdir(vdir):
        ### Use run10x default location
        vdir = os.path.join(os.path.dirname(outs_dir), "velocyto")
    if os.path.isdir(vdir):
        import scvelo as scv
        ld = [x for x in os.listdir(vdir) if x.endswith(".loom")]
        if len(ld) != 1:
            warn("number of loom in %s is not 1" % vdir)
        else:
            with sw("Reading loom"):
                vdata = scv.read_loom(os.path.join(vdir, ld[0]))
            vdata.obs.index = [x.split(":")[1] + "-1" for x in vdata.obs_names]
            vdata.var_names_make_unique()
    with sw("Reading H5"):
        adata = sc.read_10x_h5(h5, gex_only=True)
        adata.var_names_make_unique()
    if tss is not None and os.path.exists(tss):
        tss = get_tss(tss)
        adata.var["interval"] = [tss["interval"].get(g, "NA") for g in adata.var["gene_ids"]]
    if vdata is not None:
        with sw("Adding spliced/unspliced counts"):
            I = np.intersect1d(vdata.obs_names, adata.obs_names)
            adata = adata[I, vdata.var_names].copy()
            adata.layers = vdata[I, :].layers.copy()
    if os.path.isfile(bcfile):
        filtered_barcodes = pd.read_csv(bcfile, header=None, sep="\t")[0].values
        adata.obs["filtered"] = adata.obs_names.isin(filtered_barcodes).astype(str)
    if sample is not None:
        adata.obs.index = ["%s#%s" % (sample, bc) for bc in adata.obs_names]
        adata.obs[kwargs.get("sample_name", "Sample")] = sample
    with sw("Converting counts to int"):
        adata.X = adata.X.astype(np.int32)
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "Rps", "Rpl"))
    with sw("Calculating QC metrics"):
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, percent_top=[])
    with sw("Writing H5AD"):
        adata.uns = benj.convert_dict(adata.uns)
        adata.write_h5ad(output, compression="gzip", compression_opts=compression)

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--sample-name", default="Sample")
    ap.add_argument("--output", required=True)
    ap.add_argument("--compression", type=int, default=9)
    ap.add_argument("--tss", help="TSS.bed file from refdata-cellranger-arc-*/regions/tss.bed")
    ap.add_argument("--use-muon", dest="use_muon", action="store_true")
    ap.add_argument("--no-use-muon", dest="use_muon", action="store_false")
    ap.set_defaults(use_muon=True)
    args = vars(ap.parse_args())
    run(**args)
