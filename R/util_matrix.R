

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
    u = u[!is.na(u)]
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
    if (rows) {
        return(mat[ord,])
    } else {
        return(mat[,ord])
    }
}

#' Diagonally cluster the matrix
#'
#' This function clusters a matrix
#'
#' @param mat Input matrix
#' @return A list of the sorted matrix and associated columns
#' @export
diag.mat3 <- function(mat, ratio=0.5, cutoff=0.25, rows=TRUE) {
    prev.cutoff = attr(mat, "cutoff")
    if (rows) {
        ord = order(rowSums(mat > cutoff)  > ratio * ncol(mat),
                    apply(mat,1,which.max), decreasing=TRUE)
        mat = mat[ord,]
        cto = apply(mat, 1, which.max)
        idx = rowSums(mat > cutoff)  > ratio * ncol(mat)
        cto[idx] = 0
    } else {
        ord = order(colSums(mat > cutoff)  > ratio * nrow(mat),
                    apply(mat,2,which.max), decreasing=TRUE)
        mat = mat[,ord]
        cto = apply(mat, 2, which.max)
        idx = colSums(mat > cutoff)  > ratio * nrow(mat)
        cto[idx] = 0

    }
    attr(mat, "cutoff") = prev.cutoff
    if (is.null(attr(mat, "cutoff"))) {
        if (rows) {
            attr(mat, "cutoff") = list(row=cto, col=NULL)
        } else {
            attr(mat, "cutoff") = list(row=NULL, col=cto)
        }
    } else {
        if (rows) {
            attr(mat, "cutoff")$row = cto
        } else {
            attr(mat, "cutoff")$col = cto
        }
    }
    return(mat)
}

#' Pivot a dataframe using Matrix::sparseMatrix
#'
#' @param df Dataframe
#' @param row Row of dataframe
#' @param col Col of dataframe
#' @param val Value to be added
#' @param dense logical to determine whether to densify
#' @return pivoted matrix
#' @export
pivot <- function(df, row, col, val, dense=TRUE) {
    row = as.factor(df[[row]])
    col = as.factor(df[[col]])
    M = Matrix::sparseMatrix(i=as.integer(row),
                             j=as.integer(col),
                             x=df[[val]],
                             dimnames=list(levels(row), levels(col)),
                             dims=c(length(levels(row)), length(levels(col))))
    if (dense) {
        return(as.matrix(M))
    } else {
        return(M)
    }
}
