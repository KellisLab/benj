#!/usr/bin/env python3

if __name__ == "__main__":
    import benj
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", dest="h5ad", required=True)
    ap.add_argument("-o", "--output", type=str, default=None)
    ap.add_argument("-t", "--tsv")
    ap.add_argument("-b", "--batch", type=str, default=None)
    ap.add_argument("-p", "--plot", nargs="+")
    ap.add_argument("-r", "--resolution", default=1., type=float)
    ap.add_argument("-l", "--leiden", default="overall_clust")
    ap.add_argument("--prefix", default="C")
    ap.add_argument("--genome", required=True)
    ap.add_argument("-c", "--min-cells-per-peak", type=int, default=10, dest="min_n_cells_by_counts")
    ap.add_argument("-j", "--jaspar", default="JASPAR2022")
    ap.add_argument("-s", "--species", type=int, required=True)
    ap.add_argument("--cor-cutoff", type=float, default=0.8)
    ap.add_argument("--no-use-harmony", dest="use_harmony", action="store_false")
    ap.add_argument("--use-harmony", dest="use_harmony", action="store_true")
    ap.add_argument("--no-use-bbknn", dest="use_bbknn", action="store_false")
    ap.add_argument("--use-bbknn", dest="use_bbknn", action="store_true")
    ap.add_argument("--qc-cols", nargs="+", default=["log1p_total_counts"])
    ap.add_argument("--min-dist", type=float, default=0.3)
    ap.add_argument("--max-iter-harmony", type=int, default=50)
    ap.add_argument("--compression", type=int, default=6)
    ap.set_defaults(use_harmony=False, use_bbknn=True)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    adata = benj.parse_anndata(**args)
    benj.integrate_atac(adata, **args)
