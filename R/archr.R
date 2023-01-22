
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
