#' W statistics generator based on the beta(lambda) from a split LASSO path
#' @description
#' generates the split knockoff statistics W based on the beta(lambda) from
#' a split LASSO path in the intercepetion assignment step.
#'
#' @param X the design matrix
#' @param D the linear transform
#' @param y the response vector
#' @param nu the parameter for variable splitting
#' @param option options for creating the Knockoff statistics
#' option$lambda: the choice of lambda for the path
#' @return the split knockoff statistics W and various intermedia statistics
#' @export
#'
#' @importFrom glmnet glmnet
sk.W_path <- function(X, D, y, nu, option)
  # input argument:
  # X : the design matrix
  # y : the response vector
  # D : the linear transform
  # nu: the parameter for variable splitting
  # option: options for creating the Knockoff statistics
  #	option$lambda: the choice of lambda for the path
{
  m <- nrow(D)
  p <- ncol(D)
  n <- nrow(X)

  # split the dataset
  rand_rank = option$rand_rank

  ind1 = rand_rank[1:floor(n*option$frac)]
  ind2 = rand_rank[floor(n*option$frac+1):length(rand_rank)]
  X_1 = X[ind1, ]
  y_1 = y[ind1]
  X_2 = X[ind2, ]
  y_2 = y[ind2]

  ############ step 0 #############

  opts <- list(data = NA)
  opts$copy <- 'false'
  opts <- opts[-1]

  creat.result  <- sk.create(X_1, y_1, D, nu, opts)
  A_beta  <- creat.result$A_beta
  A_gamma <- creat.result$A_gamma
  tilde_y <- creat.result$tilde_y


  # set lambda

  if (is.null(option$lambda)){
    lambda_vec <- 10.^seq(0, -6, by = -0.01)
  }else{
    lambda_vec <- option$lambda
  }
  nlambda <- length(lambda_vec)

  #nlambda <- 1000
  #lambda_max = max(abs(t(cbind(A_beta,A_gamma)) %*% tilde_y)) / nrow(tilde_y)
  #lambda_min = lambda_max / 2e3
  #k = (0:(nlambda-1)) / nlambda
  #lambda_vec = lambda_max * (lambda_min/lambda_max)^k

  # set penalty
  penalty <- matrix(1,m+p,1)

  for (i in 1: p) {
    penalty[i, 1] = 0
  }

  # lasso path settings for glmnet

  fit_step0 = glmnet(cbind(A_beta,A_gamma) , tilde_y, lambda = lambda_vec, penalty.factor = penalty)
  coefs <- fit_step0$beta

  # store beta(lambda)
  betas = coefs[1:p, ]

  ############ step 1 #############

  opts <- list(data = NA)
  opts$copy <- 'true'
  opts <- opts[-1]

  creat.result  <- sk.create(X_2, y_2, D, nu, opts)
  A_beta  <- creat.result$A_beta
  A_gamma <- creat.result$A_gamma
  tilde_y <- creat.result$tilde_y
  tilde_A_gamma <- creat.result$tilde_A_gamma

  ## coef1 = coefs[(p+1):(p+m),]

  coef1 <- matrix(0,m,1)
  for (i in 1:nlambda) {
    # take beta_lambda, gamma_lambda as calculated in step 1
    y_new <- tilde_y - A_beta %*% betas[, i]
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
  if (option$m > m){
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
    # take beta_lambda, gamma_lambda as calculated in step 1
    y_new <- tilde_y - A_beta %*% betas[, i]
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
  if (option$m > m){
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

  if (is.null(option$gamma_supp)){
    return(list(W = W,
                r = r,
                Z = Z,
                t_r = t_r,
                t_Z = t_Z,
                Ws = Z * sign(Z - t_Z),
                Wst = Z * sign(Z - Z_tilde),
                Wbc = pmax(Z, t_Z) * sign(Z - t_Z),
                Wbct = pmax(Z, Z_tilde) * sign(Z - Z_tilde)))
  }else{
    return(list(gamma_supp = option$gamma_supp,
                W = W,
                r = r,
                Z = Z,
                t_r = t_r,
                t_Z = t_Z,
                Ws = Z * sign(Z - t_Z),
                Wst = Z * sign(Z - Z_tilde),
                Wbc = pmax(Z, t_Z) * sign(Z - t_Z),
                Wbct = pmax(Z, Z_tilde) * sign(Z - Z_tilde)))
  }
}
