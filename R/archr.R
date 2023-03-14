
#'
#' @export
loadAnnotationFromRNA <- function(path, sep="\t") {
    df = read.table(path, sep=sep, header=TRUE, row.names=1, comment.char="")
    spt = strsplit(rownames(df), "#")
    df$bc = sapply(spt, "[[", 1)
    df$Sample = sapply(spt, "[[", 2)
    rownames(df) = paste0(df$Sample, "#", df$bc)
    return(df)
}

#' @export
extractMarkerPeaks <- function(markerPeaks, FDR=0.10, Log2FC=0.5) {
    rD = as.data.frame(SummarizedExperiment::rowData(markerPeaks))
    A = SummarizedExperiment::assays(markerPeaks)
    mFilter = (A$FDR < FDR) & (A$Log2FC >= Log2FC)
    ms = Matrix::summary(Matrix::Matrix(mFilter))
    rD$col = NA
    rD$col[ms$i] = colnames(markerPeaks)[ms$j]
    rD = rD[!is.na(rD$col),]
    return(with(rD, GenomicRanges::GRanges(seqnames=seqnames, ranges=IRanges::IRanges(start, end), col=col)))
}

#' Create annotation from GTF or GFF file
#'
#' @param gff Path to GTF or GFF file
#' @param OrgDb org.Hs.eg.db::org.Hs.eg.db or similar object
#' @param dataSource String to create TxDb, e.g. "GENCODEv43"
#' @param organism Organism string, e.g. Homo sapiens
#' @param annoStyle How to join gene ids
#' @return ArchR style createGeneAnnotation object
#' @export
createGeneAnnotationGFF <- function(gff, OrgDb, dataSource, organism, annoStyle="ENSEMBL") {
    gdf = gff3_symbols(gff)
    gdf$orig_gene_id = gdf$gene_id
    gdf$gene_id = gsub("[.][0-9]+$", "", gdf$gene_id)
    rownames(gdf) = gdf$gene_id
    txdb = GenomicFeatures::makeTxDbFromGFF(gff, dataSource=dataSource, organism=organism)
    cga = ArchR::createGeneAnnotation(TxDb=txdb,
                                      OrgDb=OrgDb,
                                      annoStyle=annoStyle)
    return(S4Vectors::SimpleList(lapply(cga, function(obj) {
        if (all(c("gene_id", "symbol") %in% colnames(S4Vectors::mcols(obj)))) {
            S4Vectors::mcols(obj)$symbol = gdf[S4Vectors::mcols(obj)$gene_id, "gene_name"]
        }
        return(obj)
    })))
}
