

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

#' @export
gff3_symbols <- function(gff) {
    gff = rtracklayer::readGFF(gff)
    gff = gff[gff$type == "gene",]
    gff$gene_name = var_names_make_unique(gff$gene_name)
    return(gff)
}


#' Promoter enrichment of a gene set
#'
#' @param gene_list Vector of genes
#' @param universe Gene universe to compare against
#' @param upstream Base pairs upstream for promoters
#' @param downstream Base pairs downstream for promoters
#' @param jaspar JASPAR object for motifs. Default is JASPAR2020::JASPAR2020
#' @param gff GTF or GFF file used for gene locations (e.g. from Gencode)
#' @param ... Passed in benj::motif_overlap
#' @return dataframe
#' @export
promoter_enrichment <- function(gene_list, universe=NULL, upstream=2000, downstream=200,
                                jaspar=NULL,
                                gff="/home/Genomes/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz", ...) {
    gff = as.data.frame(rtracklayer::readGFF(gff))
    gff = gff[gff$type == "gene",]
    if (!is.null(universe)) {
        gff = gff[gff$gene_name %in% universe,]
    }
    gr = with(gff,
              GenomicRanges::GRanges(seqnames=seqid,
                                     ranges=IRanges::IRanges(start, end),
                                     gene_name=gene_name,
                                     gene_id=gene_id,
                                     strand=strand,
                                     gene_type=gene_type))
    names(gr) = benj::var_names_make_unique(gr$gene_name)
    prom.gr = GenomicRanges::promoters(gr, upstream=upstream, downstream=downstream)
    prom.gr$to_enrich = names(gr) %in% gene_list
    if (is.null(jaspar)) {
        jaspar = JASPAR2020::JASPAR2020
    }
    mo = motif_overlap(prom.gr, jaspar, ...)
    return(enrich_se(mo, groupby="to_enrich"))
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
