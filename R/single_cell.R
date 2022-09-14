
#' Add rowRanges to SingleCellExperiment from peak-style rownames
#'
#' Takes peak ranges in the form of chrom:start-end
#' to extract genomicranges
#'
#' @param sce SingleCellExperiment object
#' @return A SingleCellExperiment with rowRanges representing the peaks
#' @export
add_peak_row_ranges <- function(sce) {
    peak_names = rownames(sce)
    S = strsplit(peak_names, "[-:]")
    df = data.frame(chrom = sapply(S, "[[", 1),
                    begin = as.integer(sapply(S, "[[", 2)),
                    end = as.integer(sapply(S, "[[", 3)),
                    row.names=peak_names)
    rowRanges(sce) = with(df,
                          GenomicRanges::GRanges(seqnames=chrom,
                                                 ranges=IRanges::IRanges(begin, end)))
    names(rowRanges(sce)) = peak_names
    return(sce)
}
