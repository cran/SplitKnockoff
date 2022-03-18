# R function for filter of structural sparsity problem
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)

# KNOCKOFFS.PRIVATE.DECOMPOSE  Decompose design matrix X for knockoff creation

#' make SVD as well as orthogonal complements
#'
#'
#'
#' @param X the input matrix
#' @param D the linear transformer
#'
#' @return U
#' @return S
#' @return V
#' @return U_perp : orthogonal complement for U
#' @examples
#' library(mvtnorm)
#' n = 350
#' p = 100
#' m = 200
#' Sigma = matrix(0, p, p)
#' D <- matrix(0,m,p)
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)
#' decompose.result <- sk.decompose(X,D)
#' U_perp <- decompose.result$U_perp
#' @export
#'
sk.decompose <- function(X, D){
  # Check dimensions.
  m <- nrow(D)
  n <- nrow(X)
  p <- ncol(X)
  if(n < 2*p){
    print("knockoff:DimensionError: Data matrix must have n >= 2p")
  }
  # Factorize X as X = USV' (reduced SVD).
  svd.result <- canonicalSVD(X)
  S <- svd.result$S
  S <- diag(S)
  V <- svd.result$V
  U <- svd.result$U
  UU <- cbind(U,matrix(0,n,m))
  qrresult <- qr(UU)
  Qreslt <- qr.Q(qrresult)
  U_perp <- Qreslt[,(p+1):(p+m)]
  structure(list(call = match.call(),
                 U = U,
                 S = S,
                 V = V,
                 U_perp =U_perp),
            class = 'decompose.result')
}

