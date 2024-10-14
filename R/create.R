#' generate split knockoff copies
#' @description
#' Gives the variable splitting design matrix
#' and response vector. It will also create a
#' split knockoff copy if required.
#'
#' @param X the design matrix
#' @param y the response vector
#' @param D the linear transform
#' @param nu the parameter for variable splitting
#' @param option options for creating the Knockoff copy
#' option$copy true : create a knockoff copy;
#'
#' @return A_beta: the design matrix for beta after variable splitting
#' @return A_gamma: the design matrix for gamma after variable splitting
#' @return tilde_y: the response vector after variable splitting.
#' @return tilde_A_gamma: the knockoff copy of A_beta; will be NULL if option$copy = false.
#' @export
#'
#' @examples
#' option <- list()
#' option$q <- 0.2
#' option$method <- 'knockoff'
#' option$normalize <- 'true'
#' option$lambda <- 10.^seq(0, -6, by=-0.01)
#' option$nu <- 10
#' option$copy <- 'true'
#' library(mvtnorm)
#' sigma <-1
#' p <- 100
#' D <- diag(p)
#' m <- nrow(D)
#' n <- 350
#' nu = 10
#' c = 0.5
#' Sigma = matrix(0, p, p)
#' for( i in 1: p){
#'   for(j in 1: p){
#'     Sigma[i, j] <- c^(abs(i - j))
#'  }
#' }
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)
#' beta_true <- matrix(0, p, 1)
#' varepsilon <- rnorm(n) * sqrt(sigma)
#' y <- X %*% beta_true + varepsilon
#' creat.result  <- sk.create(X, y, D, nu, option)
#' A_beta  <- creat.result$A_beta
#' A_gamma <- creat.result$A_gamma
#' tilde_y <- creat.result$tilde_y
#' tilde_A_gamma <- creat.result$tilde_A_gamma
sk.create <- function(X, y, D, nu, option)
  # sk.create gives the variable splitting design matrix
  # [A_beta, A_gamma] and response vector tilde_y. It will also create a
  # split knockoff copy for A_gamma if required.
  #
  # Input Arguments:
  # X : the design matrix.
  # y : the response vector.
  # D : the linear transform.
  # nu: the parameter for variable splitting.
  # option: options for creating the Knockoff copy.
  #	option.copy = true : create a knockoff copy.
#
# Output Arguments:
# A_beta: the design matrix for beta after variable splitting.
# A_gamma: the design matrix for gamma after variable splitting.
# tilde_y: the response vector after variable splitting.
# tilde_A_gamma: the knockoff copy of A_beta; will be [] if option.copy =
# false.
{

  eps = 1e-10

  n <- nrow(X)
  m <- nrow(D)

  # calculate A_beta, A_gamma
  A_beta <- rbind(X/sqrt(n),D/sqrt(nu))
  A_gamma <- rbind(matrix(0,n,m),-diag(m)/sqrt(nu))

  # calculate tilde_y
  tilde_y <- rbind(matrix(y, ncol = 1)/sqrt(n),matrix(0,m,1))

  ##s_size <- 2-option$eta
  s_size <- 2

  if (option$copy == 'true'){
    # calculte inverse for Sigma_{beta, beta}
    Sigma_bb <- t(A_beta) %*% A_beta

    if(sum(abs(eigen(Sigma_bb)$values)>1e-6) == nrow(Sigma_bb)){
      Sigma_bb_inv = solve(Sigma_bb)
    }else{Sigma_bb_inv <- MASS::ginv(Sigma_bb)}

    # calculate Sigma_{gamma, gamma}, etc
    ##svd.result = canonicalSVD(A_gamma)
    ##S <- svd.result$S
    # S <- diag(S)
    ##U <- svd.result$U
    ##V <- svd.result$V
    # Sigma_gg = (V %*% Matrix(S^2)) %*% t(V)
    ##Sigma_gg = (V %*% diag(S^2)) %*% t(V)
    Sigma_gg = t(A_gamma) %*% A_gamma
    Sigma_gb = t(A_gamma) %*% A_beta
    Sigma_bg = t(Sigma_gb)

    # calculate C_nu
    ##H <- diag(m) - t(U) %*% A_beta %*% Sigma_bb_inv %*% t(A_beta) %*% U
    ##H <- (H + t(H))/2
    ##H_ei <- eigen(H)
    ##C = V %*% diag(S) %*% H %*% diag(S) %*% t(V)
    ##C = (C + t(C))/2
    ##C_inv = V %*% diag(1/S) %*% H_ei$vectors %*% diag(1/H_ei$values) %*% t(H_ei$vectors) %*% diag(1/S) %*% t(V)
    C = Sigma_gg - Sigma_gb %*% Sigma_bb_inv %*% Sigma_bg
    C = (C + t(C))/2
    C_inv = solve(C)
    t <- min(s_size * min(eigen(C)$values), 1/nu)
    t <- matrix(1,nrow(C),1)* t
    t <- array(t)
    # generate s
    diag_s <- diag(t)


    # calculate K^T K = 2S-S C_nu^{-1} S
    KK = 2 * diag_s - diag_s %*% C_inv %*% diag_s
    KK = (KK + t(KK))/2
    e_K <- eigen(KK)
    ##KK_lam = eigen(KK)$val
    ##KK_lam[abs(KK_lam) < eps] = 0
    ##See <- diag(KK_lam)
    See <- e_K$val
    See[See < 0] = 0
    See <- diag(See)
    Uee <- e_K$vec
    K = Uee %*% sqrt(See) %*% t(Uee)

    # calculate U=[U1;U_2] where U_2 = 0_m* m
    # U_1 is an orthogonal complement of X
    decompose.result <- sk.decompose(X, D)
    U_perp <- decompose.result$U_perp
    U_1 = U_perp[,1:m]
    U = rbind(U_1, matrix(0,m,m))

    # calculate sigma_beta beta^{-1} sigma_beta gamma
    short = Sigma_bb_inv %*% Sigma_bg


    # calculate tilde_A_gamma
    tilde_A_gamma = A_gamma %*% (diag(m) - C_inv %*% diag_s) + A_beta %*%short%*%C_inv%*%diag_s+ U %*% K
  }else{
    tilde_A_gamma = NULL
  }




  return(list(A_beta = A_beta, A_gamma = A_gamma, tilde_y = tilde_y, tilde_A_gamma = tilde_A_gamma))

}
