
#' Parse string range to seqnames, start, end
#'
#' @param string Input string in format seqnames:start-end
#' @param ranges Boolean whether to return as GRanges
#' @param name Column to keep original name as
#' @param sep Suffix to name duplicates. "-" works well as it makes idempotent
#' @return data.frame or GRanges
#' @export
parse.range <- function(string, ranges=TRUE, name=NULL, sep="-") {
    S = strsplit(string, "[:-]")
    df = data.frame(seqnames = sapply(S, "[[", 1),
                    start = as.integer(sapply(S, "[[", 2)),
                    end = as.integer(sapply(S, "[[", 3)),
                    row.names = make.unique(string, sep))
    if (ranges) {
        gr = with(df, GenomicRanges::GRanges(seqnames=seqnames,
                                             ranges=IRanges::IRanges(start, end)))
        names(gr) = rownames(df)
        if (!is.null(name)) {
            GenomicRanges::elementMetadata(gr)[[name]] = string
        }
        return(gr)
    } else {
        if (!is.null(name)) {
            df[[name]] = string
        }
        return(df)
    }
}


run.great <- function(gr, split=NULL, bg=TRUE) {
    grcol = colnames(GenomicRanges::elementMetadata(gr))
    if (!is.null(split) & all(split %in% grcol)) {

    }
}
