
def read_peaks(peaks:str, bed:bool=False, qc=[]):
    if bed:
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
    for col in qc:
        if col in pk.columns:
            pk = pd.concat([pk, pd.get_dummies(pk[col]) > 0], axis=1)
    return pk
    
def count_atac(fragments:str,
               sample:str,
               cell_metadata, peaks, output:str=None,
               compression:int=6,
               blacklist=None,
               max_value:int=127,
               stranded:bool=False,
               qc=["peakType"],
               sample_name="Sample",
               bed:bool=False,
               **kwargs):
    import os
    import pandas as pd
    import numpy as np
    import anndata
    from muon import atac as ac
    import scanpy as sc
    import pyranges
    from .timer import template as stopwatch
    from .utils import convert_dict, convert_X
    sw = stopwatch()
    if not isinstance(cell_metadata, pd.DataFrame):
        with sw("Reading cell metadata"):
            cm = pd.read_csv(cell_metadata, sep="\t", index_col=0)
    else:
        cm = cell_metadata
    cm = cm.loc[cm[sample_name] == sample, :]
    old_index = cm.index.values.astype(str)
    if kwargs.get("bc-sample", False):
        cm.index = [x[1].split("-")[0] + "-" + sample for x in cm.index.str.split("#")]
    else:
        cm.index = [x[1] for x in cm.index.str.split("#")]
    with sw("Reading peaks"):
        pk = read_peaks(peaks)
    with sw("Filtering contigs"):
        import pysam
        contigs = pysam.TabixFile(fragments, parser=pysam.asBed()).contigs
        unused = np.setdiff1d(pk["seqnames"], contigs)
        if len(unused) > 0:
            print("Warning: Unused contigs ", ",".join(unused))
        pk = pk.loc[pk["seqnames"].isin(contigs), :]
    pkbl = None
    if blacklist is not None and os.path.exists(blacklist):
        with sw("Reading blacklist"):
            bl = pd.read_csv(blacklist, sep="\t", header=None).iloc[:, [0,1,2]]
            blgr = pyranges.from_dict({"Chromosome": bl[0].values,
                                       "Start": bl[1].values,
                                       "End": bl[2].values})
            pkgr = pyranges.from_dict({"Chromosome": pk["seqnames"].values,
                                       "Start": pk["start"].values,
                                       "End": pk["end"].values,
                                       "name": pk.index.values})
            pkbl = pkgr.intersect(blgr).df
            badnames = pd.unique(pkbl["name"])
            if len(badnames) > 0: ### Only set if bad peaks exist!
                pk["overlaps_blacklist"] = pk.index.isin(badnames)
            if "strand" in pk.columns:
                pkbl["strand"] = pk.loc[pkbl["name"].values, "strand"].values
    adata = anndata.AnnData(obs=cm)
    ac.tl.locate_fragments(adata, os.path.abspath(fragments))
    with sw("Counting peaks"):
        adata = ac.tl.count_fragments_features(adata,
                                               pk.rename({"seqnames": "Chromosome", "start": "Start", "end": "End", "strand": "Strand"}, axis=1),
                                               extend_upstream=0, extend_downstream=0, stranded=stranded)
        adata.X.sum_duplicates()
    ac.tl.locate_fragments(adata, os.path.abspath(fragments))
    if pkbl is not None and isinstance(pkbl, pd.DataFrame) and pkbl.shape[0] > 0:
        with sw("Converting raw counts to int"):
            adata.X = convert_X(adata.X)
        with sw("Removing reads in peaks overlapping blacklist"):
            import scipy.sparse
            from .utils import index_of
            bdata = ac.tl.count_fragments_features(adata, pkbl, extend_upstream=0, extend_downstream=0, stranded=stranded)
            IC = index_of(pkbl["name"], adata.var_names)
            S = scipy.sparse.csr_matrix((np.ones(len(IC)), (np.arange(len(IC)), IC)),
                                        dtype=int, shape=(bdata.shape[1], adata.shape[1]))
            adata.X = adata.X - bdata.X.dot(S)
            del bdata, S
            adata.X.sum_duplicates()
    if max_value > 0:
        with sw("Clipping peak counts to %d" % max_value):
            adata.X.data = adata.X.data.clip(0, max_value)
    with sw("Converting counts to int"):
        from .utils import convert_X
        adata.X = convert_X(adata.X)
    adata.obs.index = old_index
    if len(set(adata.var.columns) & set(["interval", "nearestGene", "distToGeneStart", "peakType"])) == 4:
        with sw("Adding peak annotation"):
            annot = adata.var.loc[:, ["interval", "nearestGene", "distToGeneStart", "peakType"]].copy()
            annot.columns = ["peak", "gene", "distance", "peak_type"]
            ac.tl.add_peak_annotation(adata, annot)
    if adata.var_names.duplicated().sum() > 0:
        ### make sure cellranger style annotation is not duplicating peaks
        adata = adata[:, ~adata.var_names.duplicated()].copy()
    adata.uns = convert_dict(adata.uns)
    with sw("Calculating QC metrics"):
        qc_vars = []
        for col in qc:
            if col in pk.columns:
                ### Make sure column exists, and each value in column is also a dummy!
                qc_vars += list(set(pk.columns) & set(pk[col]))
        sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True, percent_top=None)
    if output is not None:
        with sw("Writing to H5AD"):
            adata.write_h5ad(output, compression="gzip", compression_opts=compression)
    return adata
