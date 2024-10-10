
#' Compute Liang-Zeger cluster robust standard errors
#' for glmGamPoi::glm_gp models given a clustering variable
#' @export
glm_gp_vcovCL <- function(model, cls, verbose=FALSE) {
  X <- model$model_matrix
  stopifnot(length(cls) == nrow(X))
  res_full <- residuals(model, type="response")
  bread <- solve(t(X) %*% X)    
  meat <- array(0, dim=c(ncol(X), ncol(X), nrow(res_full)))
  C <- length(unique(cls))
  for (u in unique(cls)) {
    flag <- cls == u
    Xres <- res_full[,flag,drop=FALSE] %*% X[flag,,drop=FALSE]
    meat <- meat + vapply(seq_len(nrow(Xres)), function(i) {
      outer(Xres[i,], Xres[i,])
    }, FUN.VALUE=matrix(0, ncol(Xres), ncol(Xres)))
  }
  remove("res_full")
  V <- vapply(seq_len(nrow(res_full)), function(i) {
    bread %*% meat[,,i] %*% bread
  }, FUN.VALUE=matrix(0, ncol(X), ncol(X)))
  adj <- (C/(C-1)) * (nrow(X) - 1)/(nrow(X) - ncol(X))
  return(V * adj)
}
