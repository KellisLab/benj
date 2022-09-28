

#' Varnames to names as in Scanpy
#'
#' Converts gene names to ignore
#' duplicates so that can be made
#' indices
#'
#' @param x features
#' @return Unique features
#' @export
var_names_make_unique <- function(x) { make.unique(x, sep="-") }

gff3_symbols <- function(gff) {
    gff = rtracklayer::readGFF(gff)
    gff = gff[gff$type == "gene",]
    gff$gene_name = var_names_make_unique(gff$gene_name)
    return(gff)
}


#' TSS finding
#'
#' Alternative: use GenomicFeatures::promoters(GenomicFeatures::genes(GenomicFeatures::makeTxDbFromGff(gff, dataSource="xxx", organism="xxx xxxx")))
#' @param gff GTF or GFF file location
#' @export
tss <- function(gff) {
    gff = gff3_symbols(gff)
    gr = with(gtf, GenomicRanges::GRanges(seqnames=seqid,
                                          ranges=IRanges::IRanges(tss, tss),
                                          gene=gene_name,
                                          gene_id=gene_id,
                                          tss=tss,
                                          strand=strand))
    names(gr) = gr$gene
    return(gr)
}
