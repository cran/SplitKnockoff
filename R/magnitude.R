# R function using the method of magnitude
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)


#' split_knockoffs.statistics.magnitude.cv_mag
#'
#' calculate the CV optimal lambda in the problem
#' 1/n |y - X beta|^2 + 1/nu |D beta - gamma|^2 + lambda |gamma|_1
#' for a fixed nu.
#'
#' @param X the design matrix
#' @param y the response vector
#' @param D the linear transform
#' @param nu the parameter for variable splitting
#' @param option options for calculating cv
#' option$lambda the choice of lambda for the path
#'
#' @return lambda: CV optimal lambda
#' @usage split_knockoffs.statistics.magnitude.cv_mag(X, y, D, nu, option)
#' @export
#'
split_knockoffs.statistics.magnitude.cv_mag <- function(X, y, D, nu, option)
# input argument
# X : the design matrix
# y : the response vector
# D : the linear transform
# nu: the parameter for variable splitting
# option: options for calculating cv
# option$lambda: the choice of lambda for the path

# output argument
# lambda: CV optimal lambda
{
  k_fold = 10
  n = nrow(X)
  m = nrow(m)
  p = ncol(D)

  if (is.null(option$lambda) == 'false'){
    lambda_s = option$lambda
    nlambda = length(lambda_s)
  } else {
    power = seq(0,-10,by=-0.1)
    lambda_s = 10.^power
    nlambda = length(lambda_s)}
  test <- floor(n / k_fold)
  # randomly split data
  RNGkind("Mersenne-Twister")
  rand_rank <- sample(n)

  # create matrix to store result for split lasso test loss
  loss_sl = matrix(0,k_fold, nlambda)

  #  training set and test set
  X_train <- X
  X_train <- X_train[-test, ]
  y_train <- y
  y_train <- y_train[-test, ]
  X_test <- X[test, ]
  y_test <- y[test, ]

  # normalization
  X_train <- X_train - colMeans(X_train)
  y_train <- y_train - colMeans(y_train)
  X_test  <- X_test  - colMeans(X_test)
  y_test  <- y_test  - colMeans(y_test)

  # generate split LASSO design matrix
  n_train = nrow(X_train)
  A_beta = c(X_train/sqrt(n_train),D/sqrt(nu))
  A_gamma = c(matrix(0,n_train,m),-diag(m)/sqrt(nu))
  tilde_y = c(y_train/sqrt(n_train),matrix(0,m,1))


  # set penalty for glmnet
  penalty = matrix(0,m+p,1)
  for (tmp in 1:p) {
    penalty[tmp] = 0
  }

  # opts$standardize = false # cancel standardize for glmnet
  # opts$lambda = lambda_s # choose lambda
  # opts$intr = false
  # opts$penalty_factor = penalty # choose penalty

  # fit glmnet
  fit <- glmnet(c(A_beta, A_gamma), tilde_y,  lambda =lambda_s, penalty.factor = penalty)
  # extract beta
  coefs = fit$beta
  k=1
  for (j in 1: length(lambda_s)) {
    coef = coefs[, j]
    beta = coef[1:p]
    # calculate loss
    y_sl = X_test %*% beta
    loss_sl[k, j] = norm(y_sl-y_test)^2/test
  }
  mean_loss_sl <- colMeans(loss_sl)

  # find minimal
  lambda_number = which(mean_loss_sl == min(min(mean_loss_sl)))
  lambda = lambda_s[lambda_number]
  structure(list(call = match.call(),
                 lambda = lambda),
            class = 'lambda_result')
}


#' split_knockoffs.statistics.magnitude.W_mag
#'
#' generate the knockoff statistics W, using the method of magnitude.
#' lambda here is chosen by cross validation.
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
#' @usage split_knockoffs.statistics.magnitude.W_mag(X, D, y, nu, option)
#' @export
#'
split_knockoffs.statistics.magnitude.W_mag <- function(X, D, y, nu, option)
  # input argument:
  # X : the design matrix
  # y : the response vector
  # D : the linear transform
  # nu: the parameter for variable splitting
  # option: options for creating the Knockoff statistics
  # 	option$eta : the choice of eta for creating the knockoff copy
  # 	option$lambda: the choice of lambda for the path
  #
  # output argument
  # W: the knockoff statistics
  # Z: feature significance
  # t_Z: knockoff significance
{
  m=nrow(D)
  p=ncol(D)
  # generate the design matrix
  creat.result  <- split_knockoffs.create(X, y, D, nu, option)
  A_beta  <- creat.result$A_beta
  A_gamma <- creat.result$A_gamma
  tilde_y <- creat.result$tilde_y
  tilde_A_gamma <- creat.result$tilde_A_gamma
  ########### step 0 #############

  lambda <- split_knockoffs.statistics.magnitude.cv_mag(X, y, D, nu, option)

  fit_step0 <- glmnet(c(A_beta, A_gamma), tilde_y, lambda = lambda)
  coef <- fit_step0$beta
  beta_hat <- coef[1:p]
  y_new = tilde_y - A_beta %*% beta_hat


  ########### step 1 #############
  # opts = struct;
  # opts$lambda = lambda;
   fit_step1 <- glmnet(A_gamma, y_new, lambda = lambda)
   coef1 <- fit_step1$beta

  # calculate r and Z
  Z = abs(coef1)
  r = sign(coef1)


  ########### step 2 ############
  # lasso path settings for glmnet
  # opts = struct
  # opts$lambda = lambda

  fit_step2 <- glmnet(tilde_A_gamma, y_new, lambda = lambda)
  coef2 = fit_step2$beta

  # calculate tilde_Z tilde_r
  t_Z <- matrix(0,m, 1)
  Z_prime <- abs(coef2)
  t_r <- sign(coef2)

  for (i in 1: m) {
    if(t_r(i) == r(i)){
      # store tilde_Z when it has the same sign
      t_Z[i] = Z_prime[i]
    }
  }

  ########## W #######
  W <- max(Z, t_Z) %*% sign(Z - t_Z)
  structure(list(call = match.call(),
                 W = W,
                 Z = Z,
                 t_Z = t_Z),
            class = 'mag_result')

}
