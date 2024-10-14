#' make SVD as well as orthogonal complements
#'
#' @param A the input matrix
#' @param D the linear transform
#'
#' @return U
#' @return S
#' @return V
#' @return U_perp: orthogonal complement for U
#' @export
#'
#' @examples
#' library(mvtnorm)
#' n = 350
#' p = 100
#' D <- diag(p)
#' Sigma = matrix(0, p, p)
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)
#' decompose.result <- sk.decompose(X, D)
#' U_perp <- decompose.result$U_perp
sk.decompose <- function(A, D){
  n <- nrow(A)
  p <- ncol(A)
  m <- nrow(D)
  # Factorize A as A = USV' (reduced SVD).
  svd.result <- canonicalSVD(A)
  S <- svd.result$S
  S <- diag(S)
  V <- svd.result$V
  U <- svd.result$U
  UU <- cbind(U,matrix(0,n,m))
  qrresult <- qr(UU)
  Qreslt <- qr.Q(qrresult)
  U_perp <- Qreslt[,(p+1):(m+p)]
  return(list(U = U, S = S, V = V, U_perp = U_perp))
}
