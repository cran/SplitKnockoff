#' default normalization function for matrix
#' @description
#' normalize columns of a matrix.
#'
#' @param X the input martix
#'
#' @return Y: the output matrix
#' @export
#'
#' @examples
#' library(mvtnorm)
#' n = 350
#' p = 100
#' Sigma = matrix(0, p, p)
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)
#' X <- normc(X)
normc <- function(X){

  # Normalize columns of a matrix.
  if (is.null(dim(X)))
    X <- matrix(X, ncol = 1)

  n = nrow(X)
  p = ncol(X)
  X  <- X  - matrix(rep(colMeans(X), n), nrow = n, ncol = p, byrow = TRUE)
  factors <- 1 / sqrt(colSums(X^2))
  factors <- matrix(rep(factors,n), nrow = n, ncol = p, byrow = TRUE)
  # Y = X %*% factors(matrix(1,1,n),)    # This is used in Knockoffs.
  Y = X * factors * sqrt(n-1)
  return(Y)
}
