

#' Clip a vector/matrix by values
#'
#' @param mat a matrix
#' @param min minimum value
#' @param max maximum value
#' @export
vclip <- function(mat, min=-Inf, max=Inf) {
    mat[mat < min] = min
    mat[mat > max] = max
    return(mat)
}

#' Clip a vector/matrix by quantile
#'
#' @param mat a matrix
#' @param percentile A percentile in range [0, 1]
#' @param high Whether to clip upper range
#' @param low Whether to clip lower range
#' @return A clipped matrix
#' @export
qclip <- function(mat, percentile=0.01, high=TRUE, low=TRUE) {
    if ((percentile > 1) | (percentile < 0)) {
        warning(paste0("Percentile ", percentile, "is out of range (0, 1)"))
    }
    percentile = min(1-percentile, percentile)
    Q = quantile(mat, c(percentile, 1 - percentile), na.rm=TRUE)
    if (low) {
        mat[mat < Q[1]] = Q[1]
    }
    if (high) {
        mat[mat > Q[2]] = Q[2]
    }
    return(mat)
}

#' Collapse a matrix into pseudobulk
#'
#' @param x a matrix
#' @param u unique levels to set
#' @return collapsed pseudobulk matrix
#' @export
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

#' Reorder a matrix based on TSP solver
#'
#' @param mat a matrix
#' @param rows Whether to cluster rows or cols
#' @param method Method passed to base::dist
#' @return A matrix sorted by rows
#' @export
order.tsp <- function(mat, rows=TRUE, method="euclidean") {
    if (rows) {
        tsp = seriation::seriate(dist(mat, method=method))
    } else {
        tsp = seriation::seriate(dist(t(mat), method=method))
    }
    ord = seriation::get_order(tsp, 1)
    return(ifelse(rows, mat[ord,], mat[,ord]))
}

#' Diagonally cluster the matrix
#'
#' This function clusters a matrix
#'
#' @param mat Input matrix
#' @return A list of the sorted matrix and associated columns
#' @export
diag.mat2 <- function(mat, ratio=0.5, cutoff=0.25) {
    ord = order(colSums(mat > cutoff)  > ratio * nrow(mat),
                apply(mat,2,which.max), decreasing=TRUE)
    mat = mat[,ord]
    cto = apply(mat,2,which.max)
    idx = colSums(mat > cutoff)  > ratio * nrow(mat)
    cto[idx] = 0
    return(list(mat, colnames(mat), cto))
}
