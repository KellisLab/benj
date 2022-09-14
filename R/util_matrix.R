

qnorm <- function(mat, percentile=0.01) {
    percentile = max(0, min(1, percentile))
    percentile = min(1-percentile, percentile)
    Q = quantile(mat, c(percentile, 1 - percentile), na.rm=TRUE)
    mat[mat < Q[1]] = Q[1]
    mat[mat > Q[2]] = Q[2]
    return(mat)
}

make_pseudobulk <- function(x, u=NULL) {
    if (is.factor(x)) {
        u = levels(x)
    } else if (is.null(u)) {
        u = unique(x)
    }
    names(u) = u
    return(sapply(u, function(y) {
        1 * (x == y)
    }))
}
