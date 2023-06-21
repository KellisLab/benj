#!/usr/bin/env python3

def run(adata, output:str, genome:str, release:str="JASPAR2022", chunk_size:int=10000, compression:int=9, **kwargs):
    import pychromvar as pc
    with sw("Adding peak sequence"):
        pc.add_peak_seq(adata, genome_file=genome, delimiter=":|-")
    with sw("Adding GC bias"):
        pc.add_gc_bias(adata)
    with sw("Add background peaks"):
        pc.get_bg_peaks(adata)
    with sw("Fetching motifs"):
        from pyjaspar import jaspardb
        jargs = {"collection": "CORE", "tax_group": ["vertebrates"]}
        if "species" in kwargs:
            jargs["species"] = kwargs["species"]
        jdb_obj = jaspardb(release=release)
        motifs = jdb_obj.fetch_motifs(**jargs)
    with sw("Matching motifs"):
        pc.match_motif(adata, motifs=motifs)
    with sw("Computing deviations"):
        dev = pc.compute_deviations(adata, chunk_size=chunk_size)
    with sw("Copying information"):
        tbl = {"%s.%s" % (m.matrix_id, m.name): m for m in motifs}
        dev.obsm = adata.obsm
        dev.obsp = adata.obsp
        dev.obs = adata.obs
        dev.var["matrix_id"] = [tbl[k].matrix_id for k in dev.var_names.values]
        dev.var.index = [tbl[k].name for k in dev.var_names.values]
        dev.var_names_make_unique()
    with sw("Writing to disk"):
        dev.write_h5ad(output, compression="gzip", compression_opts=compression)

if __name__ == "__main__":
    import argparse
    import benj
    ap = argparse.ArgumentParser()
    ap.add_argument("-g", "--genome", required=True)
    ap.add_argument("-j", "--jaspar", dest="release", default="JASPAR2022")
    ap.add_argument("-s", "--species", type=int, required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-i", "--input", dest="h5ad", required=True)
    ap.add_argument("-c", "--chunk-size", type=int, default=100000)
    args = benj.parse_args(ap, ["log", "scanpy", "anndata"])
    sw = benj.stopwatch()
    with sw("Reading H5AD"):
        adata = benj.parse_anndata(**args)
    run(adata, **args)
