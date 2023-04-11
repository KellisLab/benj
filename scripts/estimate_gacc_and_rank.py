#!/usr/bin/env python3

def estimate_and_rank(adata, gtf:str,
                      min_upstream:int=1000,
                      max_upstream:int=100000,
                      min_downstream:int=1000,
                      max_downstream:int=100000,
                      gene_upstream:int=5000,
                      gene_downstream:int=0,
                      target_sum:int=1000,
                      gene_scale_factor:float=5.,
                      layer:str=None,
                      groupby:str=None,
                      method:str="wilcoxon",
                      **kwargs):
    import scanpy as sc
    from benj.gene_estimation import estimate_genes_archr
    gdata = estimate_genes_archr(adata, gtf=gtf,
                                 min_upstream=min_upstream, max_upstream=max_upstream,
                                 min_downstream=min_downstream, max_downstream=max_downstream,
                                 gene_upstream=gene_upstream, gene_downstream=gene_downstream,
                                 target_sum=target_sum, gene_scale_factor=gene_scale_factor,
                                 layer=layer)
    sc.pp.log1p(gdata)
    if groupby is not None:
        sc.tl.rank_genes_groups(gdata, groupby=groupby, method=method, pts=True)
    return gdata

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", dest="h5ad", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--groupby", default=None)
    ap.add_argument("--max-upstream", type=int, default=100000)
    ap.add_argument("--max-downstream", type=int, default=100000)
    ap.add_argument("--min-upstream", type=int, default=1000)
    ap.add_argument("--min-downstream", type=int, default=1000)
    ap.add_argument("--gene-upstream", type=int, default=5000)
    ap.add_argument("--gene-downstream", type=int, default=0)
    ap.add_argument("--target-sum", type=int, default=10000)
    ap.add_argument("--gene-scale-factor", type=float, default=5.)
    ap.add_argument("--method", default="wilcoxon")
    ap.add_argument("--layer", type=str, default=None)
    ap.add_argument("--compression", type=int, default=9)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    adata = benj.parse_anndata(**args)
    gdata = estimate_and_rank(adata, **args)
    gdata.write_h5ad(args["output"], compression="gzip", compression_opts=args.get("compression", 9))
