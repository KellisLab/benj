#!/usr/bin/env python3

if __name__ == "__main__":
    import benj
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", dest="h5ad", required=True)
    ap.add_argument("-o", "--output")
    ap.add_argument("-t", "--tsv")
    ap.add_argument("-b", "--batch", type=str, default=None)
    ap.add_argument("-p", "--plot", nargs="+")
    ap.add_argument("-r", "--resolution", default=1., type=float)
    ap.add_argument("-l", "--leiden", default="overall_clust")
    ap.add_argument("--prefix", default="C")
    ap.add_argument("--hvg", default=0, type=int)
    ap.add_argument("--no-use-combat", dest="use_combat", action="store_false")
    ap.add_argument("--use-combat", dest="use_combat", action="store_true")
    ap.add_argument("--no-use-harmony", dest="use_harmony", action="store_false")
    ap.add_argument("--use-harmony", dest="use_harmony", action="store_true")
    ap.add_argument("--no-use-bbknn", dest="use_bbknn", action="store_false")
    ap.add_argument("--use-bbknn", dest="use_bbknn", action="store_true")
    ap.add_argument("--no-use-scaling", dest="use_scaling", action="store_false")
    ap.add_argument("--use-scaling", dest="use_scaling", action="store_true")
    ap.add_argument("--dotplot", nargs="+")
    ap.add_argument("--target-sum", type=int, default=0)
    ap.add_argument("--celltypist")
    ap.add_argument("--compression", type=int, default=6)
    ap.add_argument("--min-dist", type=float, default=0.3)
    ap.add_argument("--no-rgg", "--no-rank-genes", dest="use_rgg", action="store_false")
    ap.add_argument("--rgg", "--rank-genes", dest="use_rgg", action="store_true")
    ap.add_argument("--rgg-tsv", type=str, default=None)
    ap.add_argument("--max-iter-harmony", type=int, default=50)
    ap.add_argument("--leiden-n-iterations", type=int, default=-1)
    ap.add_argument("--save-data", dest="save_data", action="store_true")
    ap.add_argument("--no-save-data", dest="save_data", action="store_false")
    ap.set_defaults(use_combat=False, use_harmony=False, use_bbknn=True, use_scaling=False, use_rgg=True, save_data=False)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    print(args)
    adata = benj.parse_anndata(**args)
    benj.integrate_rna(adata, **args)
