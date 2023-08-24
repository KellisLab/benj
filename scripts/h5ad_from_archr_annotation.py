#!/usr/bin/env python3

def run(fragments, sample, cell_metadata, peaks, output, compression:int=9, blacklist=None, max_value:int=127, **kwargs):
    import os
    import pandas as pd
    import numpy as np
    import anndata
    from muon import atac as ac
    import scanpy as sc
    import pyranges
    import benj
    sw = benj.stopwatch()
    with sw("Reading cell metadata"):
        cm = pd.read_csv(cell_metadata, sep="\t", index_col=0)
    cm = cm.loc[cm["Sample"] == sample, :]
    old_index = cm.index.values.astype(str)
    if kwargs.get("bc-sample", False):
        cm.index = [x[1].split("-")[0] + "-" + sample for x in cm.index.str.split("#")]
    else:
        cm.index = [x[1] for x in cm.index.str.split("#")]
    with sw("Reading peaks"):
        if kwargs.get("bed", False):
            pk = pd.read_csv(peaks, sep="\t", header=None).iloc[:, [0,1,2]]
            pk.columns = ["seqnames", "start", "end"]
        else:
            pk = pd.read_csv(peaks, sep="\t")
        if pk.shape[1] == 6 and "peak_type" in pk.columns:
            ### cellranger peak_annotation.tsv format
            pk.columns = ["seqnames", "start", "end", "nearestGene", "distToGeneStart", "peakType"]
        if "distToTSS" in pk.columns and "distToGeneStart" not in pk.columns:
            pk = pk.rename({"distToTSS": "distToGeneStart"}, axis=1)
        pk.index = ["%s:%d-%d" % (chrom, begin, end) for chrom, begin, end in zip(pk["seqnames"], pk["start"], pk["end"])]
        pk["interval"] = pk.index.values.astype(str)
        if "peakType" in pk.columns:
            pk = pd.concat([pk, pd.get_dummies(pk["peakType"]) > 0], axis=1)
    ### TODO filter out bad
    with sw("Filtering contigs"):
        import pysam
        contigs = pysam.TabixFile(fragments, parser=pysam.asBed()).contigs
        unused = np.setdiff1d(pk["seqnames"], contigs)
        if len(unused) > 0:
            print("Warning: Unused contigs ", ",".join(unused))
        pk = pk.loc[pk["seqnames"].isin(contigs), :]
    if blacklist is not None and os.path.exists(blacklist):
        with sw("Reading blacklist"):
            bl = pd.read_csv(blacklist, sep="\t", header=None).iloc[:, [0,1,2]]
            bl = pyranges.from_dict({"Chromosome": bl[0].values,
                                     "Start": bl[1].values,
                                     "End": bl[2].values})
            pkgr = pyranges.from_dict({"Chromosome": pk["seqnames"].values,
                                       "Start": pk["start"].values,
                                       "End": pk["end"].values,
                                       "name": pk.index.values})
            badnames = pd.unique(pkgr.join(bl).df["name"])
            if len(badnames) > 0: ### Only set if bad peaks exist!
                pk["blacklist"] = pk.index.isin(badnames)
    adata = anndata.AnnData(obs=cm)
    ac.tl.locate_fragments(adata, fragments)
    with sw("Counting peaks"):
        adata = ac.tl.count_fragments_features(adata, pk.rename({"seqnames": "Chromosome", "start": "Start", "end": "End"}, axis=1), extend_upstream=0, extend_downstream=0)
    if max_value > 0:
        with sw("Clipping peak counts to %d" % max_value):
            adata.X.data = adata.X.data.clip(0, max_value)
    ac.tl.locate_fragments(adata, fragments)
    with sw("Converting counts to int"):
        adata.X = benj.convert_X(adata.X)
    adata.obs.index = old_index
    if len(set(adata.var.columns) & set(["interval", "nearestGene", "distToGeneStart", "peakType"])) == 4:
        with sw("Adding peak annotation"):
            annot = adata.var.loc[:, ["interval", "nearestGene", "distToGeneStart", "peakType"]].copy()
            annot.columns = ["peak", "gene", "distance", "peak_type"]
            ac.tl.add_peak_annotation(adata, annot)
    if adata.var_names.duplicated().sum() > 0:
        ### make sure cellranger style annotation is not duplicating peaks
        adata = adata[:, ~adata.var_names.duplicated()].copy()
    adata.uns = benj.convert_dict(adata.uns)
    with sw("Calculating QC metrics"):
        qc_vars = []
        if "peakType" in pk.columns:
            qc_vars = list(set(pk.columns) & set(pk["peakType"]))
        if "blacklist" in pk.columns:
            qc_vars.append("blacklist")
        sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True, percent_top=None)
    if "blacklist" in qc_vars:
        with sw("Zeroing blacklisted peaks and recalculating QC metrics"):
            qc_vars.remove("blacklist")
            I = np.ravel(np.where(adata.var["blacklist"].values))
            adata.X[:, I] = 0
            adata.X.eliminate_zeros()
            sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True, percent_top=None)
    with sw("Writing to H5AD"):
        adata.write_h5ad(output, compression="gzip", compression_opts=compression)
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
    ap.add_argument("-b", "--blacklist")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--max-value", type=int, default=127)
    ap.add_argument("--compression", type=int, default=9)
    ap.set_defaults(bed=False)
    args = vars(ap.parse_args())
    run(**args)
