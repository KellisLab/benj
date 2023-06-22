
#' @export
readLdscMatrix <- function(h2.dir, sep="[.]", collapse="_", verbose=TRUE) {
    FL = list.files(h2.dir, full.names=TRUE, pattern=paste0(".+", sep, ".+[.]results$"))
    L = do.call(rbind, lapply(FL, function(x) {
        bname = gsub("[.]results$", "", basename(x))
        ss = strsplit(bname, sep)[[1]]
        if ((length(ss) >= 2) & (length(readLines(x)) >= 1)) {
            if (verbose) {
                print(x)
            }
            df = read.table(x, sep="\t", header=TRUE)
            df$peaks = paste0(ss[1:(length(ss)-1)], collapse=collapse)
            df$gwas = ss[length(ss)]
            return(df)
        } else {
            return(data.frame())
        }
    }))
    L$mlog10p = -log10(L$Enrichment_p)
    return(L)
}
