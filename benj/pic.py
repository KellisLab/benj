
def pic_count(fragments, peak_bed, extend:int=5):
    import os
    import numpy as np
    from scipy.sparse import csr_matrix
    import pandas as pd
    import pyranges
    import anndata
    df = pd.read_csv(fragments, sep="\t", header=None, comment="#").iloc[:, range(4)]
    df.columns = ["seqnames", "start", "end", "barcode"]
    bc = pd.Index(pd.unique(df["barcode"]))
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
    S = csr_matrix((np.ones(pic.shape[0], dtype=np.int16), (pic["barcode_index"].values, pic["peak_index"].values)),
                   shape=(len(bc), pf.shape[0]), dtype=np.int16)
    S.sum_duplicates()
    uns = {"files": {"fragments": os.path.abspath(fragments)}}
    return anndata.AnnData(S, obs=pd.DataFrame(index=bc), var=pf, uns=uns)

def pic_create_h5ad(fragments, peak_bed, sample, h5ad,
                    tss_bed:str=None,
                    blacklist_bed:str=None,
                    compression:int=6,
                    sample_name="Sample",
                    extend:int=5,
                    promoter_upstream:int=2000,
                    promoter_downstream:int=100):
    import pandas as pd
    import scanpy as sc
    adata = pic_count(fragments=fragments, peak_bed=peak_bed, extend=extend)
    qc_vars = []
    if blacklist_bed is not None:
        import pyranges
        gr = pyranges.from_dict({"Chromosome": adata.var["seqnames"].values,
                                 "Start": adata.var["start"].values,
                                 "End": adata.var["end"].values,
                                 "interval": adata.var_names})
        bl = pd.read_csv(blacklist_bed, sep="\t", header=None)
        bl = pyranges.from_dict({"Chromosome": bl[0].values,
                                 "Start": bl[1].values,
                                 "End": bl[2].values})
        bl_peaks = pd.unique(gr.join(bl).df["interval"])
        adata.var["non_blacklist"] = ~adata.var_names.isin(bl_peaks)
        qc_vars.append("non_blacklist")
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
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inpalce=True, percent_top=None)
    if sample is not None:
        adata.obs[sample_name] = sample
        adata.obs_names = [f"{sample}#{bc}" for bc in adata.obs_names]
    adata.write_h5ad(h5ad, compression="gzip", compression_opts=compression)
    return adata
