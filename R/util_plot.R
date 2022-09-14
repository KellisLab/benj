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
