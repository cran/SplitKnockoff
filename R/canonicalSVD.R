#' singular value decomposition
#' @description
#' Computes a reduced SVD without sign ambiguity. Our convention is that
#' the sign of each vector in U is chosen such that the coefficient
#' with largest absolute value is positive.
#'
#' @param X the input matrix
#'
#' @return S
#' @return U
#' @return V
#' @export
#'
#' @examples
#' nu = 10
#' n = 350
#' m = 100
#' A_gamma <- rbind(matrix(0,n,m),-diag(m)/sqrt(nu))
#' svd.result = canonicalSVD(A_gamma)
#' S <- svd.result$S
#' S <- diag(S)
#' V <- svd.result$V
canonicalSVD <- function(X){

  svd<-svd(X)
  S <- svd$d
  U <- svd$u
  V <- svd$v
  for (j in 1:min(dim(X))) {
    i <- which.max(abs(U[,j]))
    if (U[i,j] < 0){
      U[,j] = -U[,j]
      V[,j] = -V[,j]
    }
  }
  return(list(S = S, U = U, V = V))
}
