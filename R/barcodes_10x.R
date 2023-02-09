
#' @export
load_multiome_barcodes <- function() {
    atac = system.file("extdata", "10x", "Single_Index_Kit_N_Set_A.csv", package="benj")
    rna = system.file("extdata", "10x", "Dual_Index_Kit_TT_Set_A.csv", package="benj")
    atac = read.csv(atac, comment.char="#", header=FALSE)
    rna = read.csv(rna, comment.char="#", header=TRUE)
}

#' @export
download_arc <- function(outdir, name="pbmc_granulocyte_sorted_10k", force=FALSE, method="wget", ...) {
    base = paste0("https://cf.10xgenomics.com/samples/cell-arc/2.0.0/", name, "/", name, "_")
    ext.files = c("filtered_feature_bc_matrix.h5", "atac_fragments.tsv.gz", "atac_fragments.tsv.gz.tbi", "atac_peak_annotation.tsv")
    sapply(ext.files, function(x) {
        if (force | !file.exists(x)) {
            download.file(paste0(base, x), paste0(outdir, "/", x), method=method, ...)
        }
    })
}
