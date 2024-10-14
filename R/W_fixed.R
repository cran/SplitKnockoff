#' W statistics generator based on a fixed beta(lambda) = hat beta
#' @description
#' generates the split knockoff statistics W based on a fixed beta(lambda) = hat beta in the
#' intercepetion assignment step.
#'
#' @param X the design matrix
#' @param D the linear transform
#' @param y the response vector
#' @param nu the parameter for variable splitting
#' @param option options for creating the Knockoff statistics
#' option$lambda: the choice of lambda for the path
#' option$beta_hat: the choice of beta(lambda) = hat beta
#' @return the split knockoff statistics W and various intermedia statistics
#' @export
#'
#' @importFrom glmnet glmnet
sk.W_fixed <- function(X, D, y, nu, option)
# sk.W_fixed generates the split knockoff
# statistics W based on a fixed beta(lambda) = hat{beta} in the
# intercepetion assignment step.
#
# input arguments:
# X : the design matrix
# y : the response vector
# D : the linear transformation
# nu: the parameter for variable splitting
# option: options for creating the Split Knockoff statistics
{
  m <- nrow(D)

  if (m == 0){
    return(list(W = rep(0, option$m),
                Z = rep(0, option$m),
                t_Z = rep(0, option$m),
                r = rep(0, option$m),
                t_r = rep(0, option$m)))
  }

  beta_hat = option$beta_hat

  ############ step 1 #############

  opts <- list(data = NA)
  opts$copy <- 'true'
  opts <- opts[-1]

  creat.result  <- sk.create(X, y, D, nu, opts)
  A_beta  <- creat.result$A_beta
  A_gamma <- creat.result$A_gamma
  tilde_y <- creat.result$tilde_y
  tilde_A_gamma <- creat.result$tilde_A_gamma

  # set lambda

  if (is.null(option$lambda)){
    lambda_vec <- 10.^seq(0, -6, by = -0.01)
  }else{
    lambda_vec <- option$lambda
  }
  nlambda <- length(lambda_vec)

  coef1 <- matrix(0,m,1)
  # take beta_lambda, gamma_lambda as calculated in step 1
  y_new <- tilde_y - A_beta %*% beta_hat
  for (i in 1:nlambda) {
    # calculate LASSO
    # opts = struct
    lambda = lambda_vec[i]
    # opts = glmnetSet(opts)
    fit_step1 <- glmnet(A_gamma, y_new, lambda=lambda)
    coef1 <- cbind(coef1,fit_step1$beta)
  }
  coef1 <- coef1[,-1]


  # calculate r and Z
  r <- rep(0,m)
  Z <- rep(0,m)

  for (i in 1: m) {
    hit <- hittingpoint(coef1[i, ], lambda_vec)
    Z[i] <- hit$Z
    r[i] <- hit$r
  }

  # deal with the case where some features have alrealdy been screened off
  if (!is.null(option$gamma_supp)){
    temp = rep(0, option$m)
    temp[option$gamma_supp] = Z
    Z = temp

    temp = rep(0, option$m)
    temp[option$gamma_supp] = r
    r = temp
  }

  ############# step 2 ############
  coef2 <- matrix(0,m,1)
  for (i in 1:nlambda) {
    # calculate LASSO
    # opts = struct
    lambda = lambda_vec[i]
    # opts = glmnetSet(opts)
    fit_step2 <- glmnet(tilde_A_gamma, y_new, lambda=lambda)
    coef2 <- cbind(coef2,fit_step2$beta)
  }
  coef2 <- coef2[,-1]
  # calculate tilde_Z tilde_r and W
  t_r <- rep(0,m)
  t_Z <- rep(0,m)

  for (i in 1: m) {
    t_hit <- hittingpoint(coef2[i, ], lambda_vec)
    t_r[i] <- t_hit$r
    t_Z[i] <- t_hit$Z
  }

  # deal with the case where some features have alrealdy been screened off
  if (!is.null(option$gamma_supp)){
    temp = rep(0, option$m)
    temp[option$gamma_supp] = t_Z
    t_Z = temp

    temp = rep(0, option$m)
    temp[option$gamma_supp] = t_r
    t_r = temp
  }

  ############ W ###########

  if (is.null(option$W))
    option$W = 'st'

  Z_tilde = t_Z * (r == t_r)

  switch(option$W,
         's' = {
           W = Z * sign(Z - t_Z)
         },
         'st' = {
           W = Z * sign(Z - Z_tilde)
         },
         'bc' = {
           W = pmax(Z, t_Z) * sign(Z - t_Z)
         },
         'bct' = {
           W = pmax(Z, Z_tilde) * sign(Z - Z_tilde)
         })

  return(list(W = W,
              r = r,
              Z = Z,
              t_r = t_r,
              t_Z = t_Z,
              Ws = Z * sign(Z - t_Z),
              Wst = Z * sign(Z - Z_tilde),
              Wbc = pmax(Z, t_Z) * sign(Z - t_Z),
              Wbct = pmax(Z, Z_tilde) * sign(Z - Z_tilde)))
}
