# R function for calculating the statistics W
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)


#' split_knockoffs.statistics.pathorder.W_fixed
#'
#' generate the knockoff statistics W for fixed beta in the intercepetion assignment step,
#' using the method of path order.
#'
#' @param X the design matrix
#' @param D the linear transform
#' @param y the response vector
#' @param nu the parameter for variable splitting
#' @param option options for creating the Knockoff statistics
#' option$eta the choice of eta for creating the knockoff copy
#' option$lambda the choice of lambda for the path
#' option$beta_choice the fixed beta for step 0
#'
#' @return W: the knockoff statistics
#' @return Z: feature significance
#' @return t_Z: knockoff significance
#' @examples
#' option <- array(data = NA, dim = length(data), dimnames = NULL)
#' option$q <- 0.2
#' option$eta <- 0.1
#' option$method <- 'knockoff'
#' option$stage0 <- 'path'
#' option$normalize <- 'true'
#' option$cv_rule <- 'min'
#' option$lambda <- 10.^seq(0, -6, by=-0.01)
#' option$nu <- 10
#' option$copy <- 'true'
#' option <- option[-1]
#' library(mvtnorm)
#' sigma <-1
#' n = 350
#' p = 100
#' nu = 10
#' Sigma = matrix(0, p, p)
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)
#' p <- 100
#' D <- diag(p)
#' beta_true <- matrix(0, p, 1)
#' varepsilon <- rnorm(n) * sqrt(sigma)
#' y <- X %*% beta_true + varepsilon
#' @export
#'
#'
split_knockoffs.statistics.pathorder.W_fixed <-function(X, D, y, nu, option)

# input argument:
# X : the design matrix
# y : the response vector
# D : the linear transform
# nu: the parameter for variable splitting
# option: options for creating the Knockoff statistics
#	option$eta : the choice of eta for creating the knockoff copy
#	option$lambda: the choice of lambda for the path
#	option$beta_choice : the fixed beta for step 0

# output argument
# W: the knockoff statistics
# Z: feature significance
# t_Z: knockoff significance
{
  m <- nrow(D)
  option$copy = 'true'

  # generate the design matrix
  result  <- split_knockoffs.create(X, y, D, nu, option)
  A_beta  <- result$A_beta
  A_gamma <- result$A_gamma
  tilde_y <- result$tilde_y
  tilde_A_gamma <- result$tilde_A_gamma

  ############ step 0 #############

  beta_hat <- option$beta_choice

  # calculate new response vector
  y_new <- tilde_y - A_beta %*% beta_hat

  # appoint a set of lambda
  lambda_vec <- option$lambda

  ############ step 1 #############
  # opts = struct
  if(is.null(lambda_vec) == 'false'){
    option$lambda <- lambda_vec
  }

  fit_step1 <- glmnet(A_gamma, y_new, standardize = FALSE)
  lambda_vec <- fit_step1$lambda
  coef1 <- fit_step1$beta

  # calculate r and Z
  r <- matrix(0, m, 1)
  Z <- matrix(0, m, 1)

  for (i in 1: m) {
    list[Z(i), r(i)] <- private.hittingpoint(coef1[i, ], lambda_vec)
  }

  ############# step 2 ############
  # lasso path settings for glmnet

  #opts = struct
  #opts.lambda = lambda_vec

  fit_step2 = glmnet(tilde_A_gamma, y_new, lambda = lambda_vec)
  coef2 <- fit_step2$beta

  # calculate tilde_Z tilde_r
  t_Z <- matrix(0,m, 1)
  t_r <- matrix(0,m, 1)

  for (i in 1: m) {
    result <- private.hittingpoint(coef2[i, ], lambda_vec)
    tilde_Z <-result$Z
    if(t_r[i] == r[i]){
      # store tilde_Z when it has the same sign
      t_Z[i] = tilde_Z
    }
  }

  ############ W ###########
  W <- max(Z, t_Z) %*% sign(Z - t_Z)
  structure(list(call = match.call(),
                 W = W,
                 Z = Z,
                 t_Z = t_Z),
            class = 'fixed_result')
}


#' split_knockoffs.statistics.pathorder.W_path
#'
#' generate the knockoff statistics W for beta from a split LASSO path
#' in the intercepetion assignment step, using the method of path order.
#'
#' @param X the design matrix
#' @param D the linear transform
#' @param y the response vector
#' @param nu the parameter for variable splitting
#' @param option options for creating the Knockoff statistics
#' option$eta the choice of eta for creating the knockoff copy
#' option$lambda the choice of lambda for the path
#'
#' @return W: the knockoff statistics
#' @return Z: feature significance
#' @return t_Z: knockoff significance
#' @examples
#' k <- 20   # sparsity level
#' A <- 1    # magnitude
#' n <- 350  # sample size
#' p <- 100  # dimension of variables
#' c <- 0.5  # feature correlation
#' sigma <-1 # noise level

#' option <- array(data = NA, dim = length(data), dimnames = NULL)
#' option$q <- 0.2
#' option$eta <- 0.1
#' option$method <- 'knockoff'
#' option$stage0 <- 'path'
#' option$normalize <- 'true'
#' option$cv_rule <- 'min'
#' option$lambda <- 10.^seq(0, -6, by=-0.01)
#' option$nu <- 10
#' option$copy <- 'true'
#' option <- option[-1]


#' # generate D
#' D <- diag(p)
#' m <- nrow(D)

#' # generate X
#' Sigma = matrix(0, p, p)
#' for( i in 1: p){
#'   for(j in 1: p){
#'     Sigma[i, j] <- c^(abs(i - j))
#'   }
#' }

#' library(mvtnorm)
#' set.seed(100)
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)

#' # generate beta and gamma
#' beta_true <- matrix(0, p, 1)
#' for( i in 1: k){
#'   beta_true[i, 1] = A
#'   if ( i%%3 == 1){
#'     beta_true[i, 1] = -A
#'   }
#' }
#' gamma_true <- D %*% beta_true

#' S0 <- which(gamma_true!=0)


#' # generate varepsilon
#' set.seed(1)

#' # generate noise and y
#' varepsilon <- rnorm(n) * sqrt(sigma)
#' y <- X %*% beta_true + varepsilon
#' nu = 10
#' @export
#'
#'
split_knockoffs.statistics.pathorder.W_path <-function(X, D, y, nu, option)
  # input argument:
  # X : the design matrix
  # y : the response vector
  # D : the linear transform
  # nu: the parameter for variable splitting
  # option: options for creating the Knockoff statistics
  #	option$eta : the choice of eta for creating the knockoff copy
  #	option$lambda: the choice of lambda for the path

  # output argument
  # W: the knockoff statistics
  # Z: feature significance
  # t_Z: knockoff significance
{
  m <- nrow(D)
  p <- ncol(D)

  # generate the design matrix
  creat.result  <- split_knockoffs.create(X, y, D, nu, option)
  A_beta  <- creat.result$A_beta
  A_gamma <- creat.result$A_gamma
  tilde_y <- creat.result$tilde_y
  tilde_A_gamma <- creat.result$tilde_A_gamma

  ############ step 0 #############

  # set lambda
  lambda_vec <- option$lambda
  nlambda <- length(lambda_vec)

  # set penalty
  penalty <- matrix(1,m+p,1)

  for (i in 1: p) {
    penalty[i, 1] = 0
  }

  # lasso path settings for glmnet

  fit_step0 = glmnet(cbind(A_beta,A_gamma) , tilde_y, lambda =lambda_vec, penalty.factor = penalty)
  coefs <- fit_step0$beta

  # store beta(lambda)
  betas = coefs[1: p, ]

  ############ step 1 #############
  coef1 = matrix(0,m, 1)
  for (i in 1: nlambda) {
    # take beta_lambda, gamma_lambda as calculated in step 1
    y_new <- tilde_y - A_beta %*% betas[, i]
    # calculate LASSO
    # opts = struct
    lambda = lambda_vec[i]
    fit_step1 = glmnet(A_gamma, y_new, lambda =lambda)
    # coef1[, i] <- fit_step1$beta
    coef <- fit_step1$beta
    coef1 <- cbind(coef1,coef)
  }
   coef1 <- coef1[,-1]
   # calculate r and Z
   r <- matrix(0,m,1)
   Z <- matrix(0,m,1)

   for (i in 1: m) {
     hit <- private.hittingpoint(coef1[i, ], lambda_vec)
     Z[i] <- hit$Z
     r[i] <- hit$r
   }



  ############# step 2 ############
   coef2 <- matrix(0,m,1)
   for (i in 1: nlambda) {
     # take beta_lambda, gamma_lambda as calculated in step 1
     y_new <- tilde_y - A_beta %*% betas[, i]
     # calculate LASSO
     # opts = struct
     lambda = lambda_vec[i]
     # opts = glmnetSet(opts)
     fit_step2 <- glmnet(tilde_A_gamma, y_new, lambda=lambda)
     coef <- fit_step2$beta
     coef2 <- cbind(coef2,coef)
   }
   coef2 <- coef2[,-1]
   # calculate tilde_Z tilde_r and W
   t_Z <- matrix(0,m, 1)
   t_r <- matrix(0,m, 1)

  for (i in 1: m) {
    t_hit <- private.hittingpoint(coef2[i, ], lambda_vec)
    t_r[i] <- t_hit$r
    if(t_r[i] == r[i]){
      # store tilde_Z when it has the same sign
      t_Z[i] = t_hit$Z
    }
  }


  ############ W ###########
   W <- matrix(0,m, 1)
   for (i in 1: m) {
   W[i] <- max(Z[i], t_Z[i]) * sign(Z[i] - t_Z[i])
   }
   structure(list(call = match.call(),
                  W = W,
                  Z = Z,
                  t_Z = t_Z),
             class = 'path_result')
}





#' split_knockoffs.statistics.sign.W_path
#'
#' a new method to generate the knockoff statistics W for beta from a split LASSO path
#' in the intercepetion assignment step, using the method of path order.
#'
#' @param X the design matrix
#' @param D the linear transform
#' @param y the response vector
#' @param nu the parameter for variable splitting
#' @param option options for creating the Knockoff statistics
#' option$eta the choice of eta for creating the knockoff copy
#' option$lambda the choice of lambda for the path
#'
#' @return W: the knockoff statistics
#' @return Z: feature significance
#' @return r: the sign estimator
#' @return t_Z: knockoff significance
#' @examples
#' k <- 20   # sparsity level
#' A <- 1    # magnitude
#' n <- 350  # sample size
#' p <- 100  # dimension of variables
#' c <- 0.5  # feature correlation
#' sigma <-1 # noise level

#' option <- array(data = NA, dim = length(data), dimnames = NULL)
#' option$q <- 0.2
#' option$eta <- 0.1
#' option$method <- 'knockoff'
#' option$stage0 <- 'path'
#' option$normalize <- 'true'
#' option$cv_rule <- 'min'
#' option$lambda <- 10.^seq(0, -6, by=-0.01)
#' option$nu <- 10
#' option$copy <- 'true'
#' option <- option[-1]


#' # generate D
#' D <- diag(p)
#' m <- nrow(D)

#' # generate X
#' Sigma = matrix(0, p, p)
#' for( i in 1: p){
#'   for(j in 1: p){
#'     Sigma[i, j] <- c^(abs(i - j))
#'   }
#' }

#' library(mvtnorm)
#' set.seed(100)
#' X <- rmvnorm(n,matrix(0, p, 1), Sigma)

#' # generate beta and gamma
#' beta_true <- matrix(0, p, 1)
#' for( i in 1: k){
#'   beta_true[i, 1] = A
#'   if ( i%%3 == 1){
#'     beta_true[i, 1] = -A
#'   }
#' }
#' gamma_true <- D %*% beta_true

#' S0 <- which(gamma_true!=0)


#' # generate varepsilon
#' set.seed(1)

#' # generate noise and y
#' varepsilon <- rnorm(n) * sqrt(sigma)
#' y <- X %*% beta_true + varepsilon
#' nu = 10
#' @export
#'
#'
split_knockoffs.statistics.sign.W_path <-function(X, D, y, nu, option)
  # input argument:
  # X : the design matrix
  # y : the response vector
  # D : the linear transform
  # nu: the parameter for variable splitting
  # option: options for creating the Knockoff statistics
  #	option$eta : the choice of eta for creating the knockoff copy
  #	option$lambda: the choice of lambda for the path

  # output argument
  # W: the knockoff statistics
# Z: feature significance
# t_Z: knockoff significance
{
  m <- nrow(D)
  p <- ncol(D)

  # generate the design matrix
  creat.result  <- split_knockoffs.create(X, y, D, nu, option)
  A_beta  <- creat.result$A_beta
  A_gamma <- creat.result$A_gamma
  tilde_y <- creat.result$tilde_y
  tilde_A_gamma <- creat.result$tilde_A_gamma

  ############ step 0 #############

  # set lambda
  lambda_vec <- option$lambda
  nlambda <- length(lambda_vec)

  # set penalty
  penalty <- matrix(1,m+p,1)

  for (i in 1: p) {
    penalty[i, 1] = 0
  }

  # lasso path settings for glmnet

  fit_step0 = glmnet(cbind(A_beta,A_gamma) , tilde_y, lambda =lambda_vec, penalty.factor = penalty)
  coefs <- fit_step0$beta

  # store beta(lambda)
  betas = coefs[1: p, ]

  ############ step 1 #############
  coef1 = coefs[(p+1):(p+m),]
  # calculate r and Z
  r <- matrix(0,m,1)
  Z <- matrix(0,m,1)

  for (i in 1: m) {
    hit <- private.hittingpoint(coef1[i, ], lambda_vec)
    Z[i] <- hit$Z
    r[i] <- hit$r
  }
  # for (i in 1: nlambda) {
  #   # take beta_lambda, gamma_lambda as calculated in step 1
  #   y_new <- tilde_y - A_beta %*% betas[, i]
  #   # calculate LASSO
  #   # opts = struct
  #   lambda = lambda_vec[i]
  #   fit_step1 = glmnet(A_gamma, y_new, lambda =lambda)
  #   # coef1[, i] <- fit_step1$beta
  #   coef <- fit_step1$beta
  #   coef1 <- cbind(coef1,coef)
  # }
  #  coef1 <- coef1[,-1]
  #
  #
  #  for (i in 1: m) {
  #    hit <- private.hittingpoint(coef1[i, ], lambda_vec)
  #    Z[i] <- hit$Z
  #    r[i] <- hit$r
  #  }
  #


  ############# step 2 ############
  coef2 <- matrix(0,m,1)
  for (i in 1: nlambda) {
    # take beta_lambda, gamma_lambda as calculated in step 1
    y_new <- tilde_y - A_beta %*% betas[, i]
    # calculate LASSO
    # opts = struct
    lambda = lambda_vec[i]
    # opts = glmnetSet(opts)
    fit_step2 <- glmnet(tilde_A_gamma, y_new, lambda=lambda)
    coef <- fit_step2$beta
    coef2 <- cbind(coef2,coef)
  }
  coef2 <- coef2[,-1]
  # calculate tilde_Z tilde_r and W
  t_Z <- matrix(0,m, 1)

  for (i in 1: m) {
    t_hit <- private.hittingpoint(coef2[i, ], lambda_vec)
    t_Z[i] <- t_hit$Z
  }


  ############ W ###########
  W <- matrix(0,m, 1)
  for (i in 1: m) {
    W[i] <- Z[i] * sign(Z[i] - t_Z[i])
  }
  structure(list(call = match.call(),
                 W = W,
                 r = r,
                 Z = Z,
                 t_Z = t_Z),
            class = 'path_result')
}

