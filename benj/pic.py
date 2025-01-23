
def pic_count(fragments, peak_bed, extend:int=5, blacklist_bed=None, max_value:int=127):
    """ PyRanges and BED are both 0-based, with [start, end) style intervals.
    """
    import os
    import numpy as np
    from scipy.sparse import csr_matrix
    import pandas as pd
    import pyranges
    import anndata
    df = pd.read_csv(fragments, sep="\t", header=None, comment="#").iloc[:, range(4)]
    df.columns = ["seqnames", "start", "end", "barcode"]
    if blacklist_bed is not None:
        gr = pyranges.from_dict({"Chromosome": df["seqnames"].values,
                                 "Start": df["start"].values,
                                 "End": df["end"].values,
                                 "frag_index": np.arange(df.shape[0])})
        bl = pd.read_csv(blacklist_bed, sep="\t", header=None, comment="#").iloc[:, range(3)]
        bl.columns = ["seqnames", "start", "end"]
        bl = pyranges.from_dict({"Chromosome": bl["seqnames"].values,
                                 "Start": bl["start"].values,
                                 "End": bl["end"].values})
        frag_index = gr.overlap(bl, invert=True).df["frag_index"].values
        del gr, bl
        df = df.iloc[frag_index, :].copy()
    bc = pd.Index(pd.unique(df["barcode"]))
    ### Subtract one for the current base pair.
    ### Add extra to LHS like GenomicRanges::resize
    lhs_extend = int((extend - 1) // 2 + (extend - 1) % 2)
    rhs_extend = int((extend - 1) // 2) + 1
    pf = pd.read_csv(peak_bed, sep="\t", header=None).iloc[:, range(3)]
    pf.columns = ["seqnames", "start", "end"]
    pf["interval"] = ["%s:%d-%d" % (seqnames, start, end) for seqnames, start, end in zip(pf.seqnames, pf.start, pf.end)]
    pf.index = pf["interval"].values
    peaks = pyranges.from_dict({"Chromosome": pf["seqnames"].values,
                                "Start": pf["start"].values,
                                "End": pf["end"].values,
                                "peak_index": np.arange(pf.shape[0])})
    left = pyranges.from_dict({"Chromosome": df["seqnames"].values,
                               "Start": df["start"].values - lhs_extend,
                               "End": df["start"].values + rhs_extend,
                               "frag_index": np.arange(df.shape[0])})
    right = pyranges.from_dict({"Chromosome": df["seqnames"].values,
                                "Start": df["end"].values - lhs_extend - 1, ### End is non-inclusive
                                "End": df["end"].values + rhs_extend - 1,
                                "frag_index": np.arange(df.shape[0])})
    ovp_s = left.join(peaks, how="left", preserve_order=True).df.loc[:, ["frag_index", "peak_index"]].sort_values(["frag_index", "peak_index"]).drop_duplicates("frag_index")
    ovp_e = right.join(peaks, how="left", preserve_order=True).df.loc[:, ["frag_index", "peak_index"]].sort_values(["frag_index", "peak_index"]).drop_duplicates("frag_index")
    del left, right, peaks
    ovp = pd.DataFrame({"frag_index": np.arange(df.shape[0]), "peak_index_s": -1, "peak_index_e": -1})
    ovp.loc[ovp_s["frag_index"].values, "peak_index_s"] = ovp_s["peak_index"].values
    ovp.loc[ovp_e["frag_index"].values, "peak_index_e"] = ovp_e["peak_index"].values
    ovp["barcode_index"] = bc.get_indexer(df["barcode"].values[ovp["frag_index"].values])
    del ovp_s, ovp_e
    pic = pd.concat((ovp.loc[ovp["peak_index_s"] != ovp["peak_index_e"], ["barcode_index", "peak_index_e"]].rename({"peak_index_e": "peak_index"}, axis=1),
                     ovp.loc[:, ["barcode_index", "peak_index_s"]].rename({"peak_index_s": "peak_index"}, axis=1)))
    pic = pic.loc[pic["peak_index"] >= 0, :]
    S = csr_matrix((np.ones(pic.shape[0], dtype=np.int64), (pic["barcode_index"].values, pic["peak_index"].values)),
                   shape=(len(bc), pf.shape[0]), dtype=np.int64)
    S.sum_duplicates()
    if max_value <= 127:
        S.data = S.data.clip(-128, max_value)
        S = S.astype(np.int8)
    elif max_value <= 32767:
        S.data = S.data.clip(-32768, max_value)
        S = S.astype(np.int16)
    elif max_value <= 2147483647:
        S.data = S.data.clip(-2147483648, max_value)
        S = S.astype(np.int32)
    uns = {"files": {"fragments": os.path.abspath(fragments)}}
    return anndata.AnnData(S, obs=pd.DataFrame(index=bc), var=pf, uns=uns)

def pic_create_h5ad(fragments, peak_bed, sample, h5ad,
                    tss_bed:str=None,
                    blacklist_bed:str=None,
                    compression:int=6,
                    sample_name="Sample",
                    extend:int=5,
                    promoter_upstream:int=2000,
                    promoter_downstream:int=100,
                    max_value=127):
    import pandas as pd
    import scanpy as sc
    adata = pic_count(fragments=fragments, peak_bed=peak_bed, extend=extend, blacklist_bed=blacklist_bed, max_value=max_value)
    qc_vars = []
    if tss_bed is not None:
        import pyranges
        from muon import atac as ac
        from .gene_estimation import get_tss
        tss = get_tss(tss_bed).rename({"left": "Start", "right": "End"}, axis=1)
        tsse = ac.tl.tss_enrichment(adata, tss, n_tss=tss.shape[0])
        gr = pyranges.from_dict({"Chromosome": adata.var["seqnames"].values,
                                 "Start": adata.var["start"].values,
                                 "End": adata.var["end"].values,
                                 "interval": adata.var_names})
        tr = pyranges.from_dict({"Chromosome": tss["Chromosome"].values,
                                 "Start": tss["Start"].values,
                                 "End": tss["End"].values})
        tss_peaks = pd.unique(gr.join(tr).df["interval"])
        adata.var["overlap_tss"] = adata.var_names.isin(tss_peaks)
        qc_vars.append("overlap_tss")
        ### Promoters
        if "strand" in tss.columns:
            tss["PromStart"] = tss["Start"] - promoter_upstream
            tss["PromEnd"] = tss["End"] + promoter_downstream
            tss.loc[tss["strand"] == "-", "PromStart"] = tss["Start"] - promoter_downstream
            tss.loc[tss["strand"] == "-", "PromEnd"] = tss["End"] + promoter_upstream
            pr = pyranges.from_dict({"Chromosome": tss["Chromosome"].values,
                                     "Start": tss["PromStart"].values,
                                     "End": tss["PromEnd"].values})
            prom_peaks = pd.unique(gr.join(pr).df["interval"])
            adata.var["overlap_promoter"] = adata.var_names.isin(prom_peaks)
            qc_vars.append("overlap_promoter")
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True, percent_top=None)
    if sample is not None:
        adata.obs[sample_name] = sample
        adata.obs_names = [f"{sample}#{bc}" for bc in adata.obs_names]
    adata.write_h5ad(h5ad, compression="gzip", compression_opts=compression)
    return adata
