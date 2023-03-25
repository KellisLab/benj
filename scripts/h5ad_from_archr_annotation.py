#!/usr/bin/env python3

def run(fragments, sample, cell_metadata, peaks, output, **kwargs):
    import pandas as pd
    import numpy as np
    import anndata
    from muon import atac as ac
    import scanpy as sc
    import benj
    sw = benj.stopwatch()
    with sw("Reading cell metadata") as _:
        cm = pd.read_csv(cell_metadata, sep="\t", index_col=0)
    cm = cm.loc[cm["Sample"] == sample, :]
    old_index = cm.index.values.astype(str)
    if kwargs.get("bc-sample", False):
        cm.index = [x[1].split("-")[0] + "-" + sample for x in cm.index.str.split("#")]
    else:
        cm.index = [x[1] for x in cm.index.str.split("#")]
    with sw("Reading peaks") as _:
        if kwargs.get("bed", False):
            pk = pd.read_csv(peaks, sep="\t", header=None).iloc[:, [0,1,2]]
            pk.columns = ["seqnames", "start", "end"]
        else:
            pk = pd.read_csv(peaks, sep="\t")
        if pk.shape[1] == 6 and "peak_type" in pk.columns:
            ### cellranger peak_annotation.tsv format
            pk.columns = ["seqnames", "start", "end", "gene", "distance", "peakType"]
        pk.index = ["%s:%d-%d" % (chrom, begin, end) for chrom, begin, end in zip(pk["seqnames"], pk["start"], pk["end"])]
        pk["interval"] = pk.index.values.astype(str)
        if "peakType" in pk.columns:
            pk = pd.concat([pk, pd.get_dummies(pk["peakType"]) > 0], axis=1)
    adata = anndata.AnnData(obs=cm)
    ac.tl.locate_fragments(adata, fragments)
    with sw("Counting peaks") as _:
        adata = ac.tl.count_fragments_features(adata, pk.rename({"seqnames": "Chromosome", "start": "Start", "end": "End"}, axis=1), extend_upstream=0, extend_downstream=0)
    ac.tl.locate_fragments(adata, fragments)
    with sw("Converting counts to int") as _:
        adata.X = adata.X.astype(np.int16)
    adata.obs.index = old_index
    if len(set(adata.var.columns) & set(["interval", "nearestGene", "distToGeneStart", "peakType"])) == 4:
        with sw("Adding peak annotation") as _:
            annot = adata.var.loc[:, ["interval", "nearestGene", "distToGeneStart", "peakType"]].copy()
            annot.columns = ["peak", "gene", "distance", "peak_type"]
            ac.tl.add_peak_annotation(adata, annot)
    adata.uns = benj.convert_dict(adata.uns)
    with sw("Calculating QC metrics") as _:
        qc_vars = []
        if "peakType" in pk.columns:
            qc_vars = list(set(pk.columns) & set(pk["peakType"]))
        sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True, percent_top=None)
    with sw("Writing to H5AD") as _:
        adata.write_h5ad(output, compression="gzip", compression_opts=9)
    return 0

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fragments", required=True)
    ap.add_argument("-s", "--sample", required=True)
    ap.add_argument("--cell-metadata", required=True)
    ap.add_argument("--peaks", required=True)
    ap.add_argument("--peaks-bed", dest="bed", action="store_true")
    ap.add_argument("--peaks-tsv", dest="bed", action="store_false")
    ap.add_argument("-o", "--output", required=True)
    ap.set_defaults(bed=False)
    args = vars(ap.parse_args())
    run(**args)
