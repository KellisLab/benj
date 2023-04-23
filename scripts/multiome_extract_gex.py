#!/usr/bin/env python3

def multiome_extract_gex(**kwargs):
    import scanpy as sc
    adata = sc.read_10x_h5(kwargs["input"], gex_only=kwargs.get("gex_only", True))
    adata.write_h5ad(kwargs["output"], compression="gzip", compression_opts=kwargs.get("compression", 9))

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--gex-only", dest="gex_only", action="store_true")
    ap.add_argument("--no-gex-only", dest="gex_only", action="store_false")
    ap.add_argument("--compression", type=int, default=9)
    ap.set_defaults(gex_only=True)
    args = vars(ap.parse_args())
    multiome_extract_gex(**args)
