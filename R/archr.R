
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
