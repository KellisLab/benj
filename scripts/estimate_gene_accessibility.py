#!/usr/bin/env python3
from typing import Iterable

def estimate_and_rank(adata, gtf:str,
                      min_upstream:int=1000,
                      max_upstream:int=100000,
                      min_downstream:int=1000,
                      max_downstream:int=100000,
                      gene_upstream:int=5000,
                      gene_downstream:int=0,
                      target_sum:int=10000,
                      gene_scale_factor:float=5.,
                      layer:str=None,
                      groupby:str=None,
                      method:str="wilcoxon",
                      plot:Iterable[str]=None,
                      celltypist:str=None,
                      over_clustering:str=None,
                      distal:bool=True,
                      log1p:bool=True,
                      tss:str=None,
                      gene:bool=True,
                      **kwargs):
    import os
    import scanpy as sc
    import pyranges
    from benj.gene_estimation import estimate_genes_archr, add_interval
    if gene:
        gdata = estimate_genes_archr(adata, gtf=gtf, feature_column=kwargs.get("feature_column"),
                                     min_upstream=min_upstream, max_upstream=max_upstream,
                                     min_downstream=min_downstream, max_downstream=max_downstream,
                                     gene_upstream=gene_upstream, gene_downstream=gene_downstream,
                                     target_sum=target_sum, gene_scale_factor=gene_scale_factor,
                                     layer=layer, log1p=log1p, distal=distal)
    else:
        from benj.count_atac import read_peaks
        feature_df = read_peaks(gtf).rename({"seqnames": "Chromosome", "start": "Start", "end": "End"}, axis=1)
        gdata = estimate_features_archr(adata, feature_df=feature_df, feature_column=kwargs.get("feature_column"),
                                        min_upstream=min_upstream, max_upstream=max_upstream,
                                        min_downstream=min_downstream, max_downstream=max_downstream,
                                        gene_upstream=gene_upstream, gene_downstream=gene_downstream,
                                        target_sum=target_sum, gene_scale_factor=gene_scale_factor,
                                        layer=layer, log1p=log1p, distal=distal)
    if tss is not None and os.path.exists(tss):
        add_interval(gdata.var, tss)
    if celltypist is not None:
        from celltypist import annotate
        if over_clustering is not None:
            ct = annotate(gdata, majority_voting=True, over_clustering=over_clustering, model=celltypist)
        else:
            ct = annotate(gdata, majority_voting=True, model=celltypist)
        for cn in ct.predicted_labels.columns:
            gdata.obs[cn] = ct.predicted_labels[cn]
    if plot is None:
        plot = []
    for item in plot:
        sc.pl.umap(gdata, color=item, save="_%s.png" % item)
    if groupby is not None:
        sc.tl.rank_genes_groups(gdata, groupby=groupby, method=method, pts=True)
    return gdata


if __name__ == "__main__":
    import argparse
    import benj
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", dest="h5ad", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--gtf", required=True, help="GTF file for gene annotation")
    ap.add_argument("--groupby", default=None, help=".obs group for ranking genes")
    ap.add_argument("--plot", nargs="+", help="Items to plot in UMAP")
    ap.add_argument("--method", default="wilcoxon", help="Method for rank_genes_groups")
    ap.add_argument("--max-upstream", type=int, default=100000)
    ap.add_argument("--max-downstream", type=int, default=100000)
    ap.add_argument("--min-upstream", type=int, default=1000)
    ap.add_argument("--min-downstream", type=int, default=1000)
    ap.add_argument("--gene-upstream", type=int, default=5000)
    ap.add_argument("--gene-downstream", type=int, default=0)
    ap.add_argument("--target-sum", type=int, default=10000)
    ap.add_argument("--gene-scale-factor", type=float, default=5.)
    ap.add_argument("--layer", type=str, default=None)
    ap.add_argument("--compression", type=int, default=6)
    ap.add_argument("--celltypist", default=None)
    ap.add_argument("--tss", help="TSS.bed file from refdata-cellranger-arc-*/regions/tss.bed used to generate \"interval\" field")
    ap.add_argument("--over-clustering", default=None, help="Overclustering column e.g. leiden used for assigning cell type")
    ap.add_argument("--log1p", dest="log1p", action="store_true")
    ap.add_argument("--no-log1p", dest="log1p", action="store_false")
    ap.add_argument("--distal", dest="distal", action="store_true")
    ap.add_argument("--no-distal", dest="distal", action="store_false")
    ap.add_argument("--gene", dest="gene", action="store_true")
    ap.add_argument("--arbitrary-feature", dest="gene", action="store_false")
    ap.add_argument("--feature-column", default="gene_id")
    ap.set_defaults(distal=True, log1p=True, gene=True)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    sw = benj.stopwatch()
    with sw("Reading H5AD"):
        adata = benj.parse_anndata(**args)
    gdata = estimate_and_rank(adata, **args)
    with sw("Writing H5AD"):
        gdata.write_h5ad(args["output"], compression="gzip", compression_opts=args.get("compression", 9))
