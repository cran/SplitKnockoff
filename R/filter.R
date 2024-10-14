#' split Knockoff filter for structural sparsity problem
#'
#' @param X the design matrix
#' @param D the response vector
#' @param y the linear transformation
#' @param option options for creating the Split Knockoff statistics.
#' option$q: the desired FDR control target.
#' option$beta: choices on beta(lambda), can be: 'path', beta(lambda)
#'       is taken from a regularization path; 'cv_beta', beta(lambda) is taken
#'       as the cross validation optimal estimator hat beta; or 'cv_all',
#'       beta(lambda) as well as nu are taken from the cross validation
#'       optimal estimators hat beta and hat nu.The default setting is
#'       'cv_all'.
#' option$lambda_cv: a set of lambda appointed for cross validation in
#'       estimating hat beta, default 10.^seq(0, -8, by = -0.4).
#' option$nu_cv: a set of nu appointed for cross validation in
#'       estimating hat beta and hat nu, default 10.^seq(0, 2, by = 0.4).
#' option$nu: a set of nu used in option.beta = 'path' or 'cv_beta' for
#'       Split Knockoffs, default 10.^seq(0, 2, by = 0.2).
#' option$lambda: a set of lambda appointed for Split LASSO path
#'       calculation, default 10.^seq(0, -6, by = -0.01).
#' option$normalize: whether to normalize the data, default true.
#' option$W: the W statistics used for Split Knockoffs, can be 's', 'st',
#' 'bc', 'bct', default 'st'.
#'
#' @return various intermedia statistics
#' @export
#'
sk.filter <- function(X, D, y, option){
  # Split Knockoff filter for structural sparsity problem.
  #
  # Input Arguments
  # X : the design matrix.
  # y : the response vector.
  # D : the linear transformation.
  # option: options for creating the Split Knockoff statistics.
  #	option.q: the desired FDR control target.
  #   option.beta: choices on \beta(\lambda), can be: 'path', \beta(\lambda)
  #       is taken from a regularization path; 'cv_beta', \beta(\lambda) is taken
  #       as the cross validation optimal estimator \hat\beta; or 'cv_all',
  #       \beta(\lambda) as well as \nu are taken from the cross validation
  #       optimal estimators \hat\beta and \hat\nu.The default setting is
  #       'cv_all'.
  #	option.lambda_cv: a set of lambda appointed for cross validation in
  #       estimating \hat\beta, default 10.*[0:-0.4:-8].
  #	option.nu_cv: a set of nu appointed for cross validation in
  #       estimating \hat\beta and \hat\nu, default 10.*[0:0.4:2].
  #   option.nu: a set of nu used in option.beta = 'path' or 'cv_beta' for
  #       Split Knockoffs, default 10.*[0:0.2:2].
  #	option.lambda: a set of lambda appointed for Split LASSO path
  #       calculation, default 10.*[0:-0.01:-6].
  #	option.normalize: whether to normalize the data, default true.
  #   option.W: the W statistics used for Split Knockoffs, can be 's', 'st',
  #       'bc', 'bct', default 'st'.

  if (option$normalize == 'true' || is.null(option$normalize)){
    X = normc(X)  # normalize(X)
    y = normc(y)  # normalize(y)
    y = y[ ,1]
  }

  X = as.matrix(X)

  q = option$q
  n = nrow(X)
  p = ncol(X)
  m = nrow(D)

  if (is.null(option$beta)){
    option$beta = 'cv_all'
  }

  if (!is.null(option$frac)){
    option$n_2 = n - floor(n * option$frac)
  }

  switch(option$beta,
         'cv_all' = {
           if (!is.null(option$seed)){
             set.seed(option$seed)
           }
           rand_rank = sample(n)

           # record random order for function W_path
           option$rand_rank = rand_rank
           ind1 = rand_rank[1: floor(n * option$frac)]
           ind2 = rand_rank[floor(n*option$frac+1):length(rand_rank)]
           X_1 = X[ind1, ]
           y_1 = y[ind1]
           X_2 = X[ind2, ]
           y_2 = y[ind2]
           n_1 = nrow(X_1)

           if (n_1 < p || (n-n_1) < (p+m)){
             # give support set estimation with dataset 1 in high
             # dimensional cases
             option$n_2 = n-n_1
             result_cv = cv_screen(X_1, y_1, D, option)
             beta_hat = result_cv$beta_hat
             stat_cv = result_cv$stat_cv
             X_new = X_2[, stat_cv$beta_supp]
             D_new = D[stat_cv$gamma_supp, stat_cv$beta_supp]
             option$gamma_supp = stat_cv$gamma_supp
           }else{
             result_cv = cv_all(X_1, y_1, D, option)
             beta_hat = result_cv$beta_hat
             stat_cv = result_cv$stat_cv
             X_new = X_2
             D_new = D
           }
           option$m = m
           option$beta_hat = beta_hat
           nu = stat_cv$nu
           stats = sk.W_fixed(X_new, D_new, y_2, nu, option)
           W = stats$W
           results = list()
           method = 'knockoff'
           results$sk = select(W, q, method)
           method = 'knockoff+'
           results$sk_plus = select(W, q, method)
           stats$nu = stat_cv$nu
           if (n_1 < p || (n-n_1) < (p+m)){
             stats$gamma_supp = stat_cv$gamma_supp
           }
         },
         'path' = {
           if (!is.null(option$seed)){
             set.seed(option$seed)
           }
           rand_rank = sample(n)

           # record random order for function W_path
           option$rand_rank = rand_rank
           ind1 = rand_rank[1: floor(n * option$frac)]
           X_1 = X[ind1, ]
           y_1 = y[ind1]
           n_1 = nrow(X_1)
           option$m = m
           if (n_1 < p || (n-n_1) < (p+m)){
             # give support set estimation with dataset 1 in high
             # dimensional cases
             option$n_2 = n-n_1
             stat_cv = cv_screen(X_1, y_1, D, option)$stat_cv
             X_new = X[, stat_cv$beta_supp]
             D_new = D[stat_cv$gamma_supp, stat_cv$beta_supp]
             option$gamma_supp = stat_cv$gamma_supp
           }else{
             X_new = X
             D_new = D
           }

           if (!is.null(option$nu)){
             nu_s = option$nu
           }else{
             nu_s = 10.^seq(0, 2, by = 0.2)
           }
           num_nu = length(nu_s)
           results = list()
           stats = list()

           for (i in 1:num_nu){
             nu = nu_s[i]
             stats[[i]] = sk.W_path(X_new, D_new, y, nu, option)
             W = stats[[i]]$W
             results[[i]] = list()
             method = 'knockoff'
             results[[i]]$sk = select(W, q, method)
             method = 'knockoff+'
             results[[i]]$sk_plus = select(W, q, method)
           }
         },
         'cv_beta' = {
           if (!is.null(option$seed)){
             set.seed(option$seed)
           }
           rand_rank = sample(n)

           # record random order for function W_path
           option$rand_rank = rand_rank
           ind1 = rand_rank[1: floor(n * option$frac)]
           ind2 = rand_rank[floor(n*option$frac+1):length(rand_rank)]
           X_1 = X[ind1, ]
           y_1 = y[ind1]
           X_2 = X[ind2, ]
           y_2 = y[ind2]
           n_1 = nrow(X_1)
           if (n_1 < p || (n-n_1) < (p+m)){
             # give support set estimation with dataset 1 in high
             # dimensional cases
             option$n_2 = n-n_1
             result_cv = cv_screen(X_1, y_1, D, option)
             beta_hat = result_cv$beta_hat
             stat_cv = result_cv$stat_cv
             X_new = X_2[, stat_cv$beta_supp]
             D_new = D[stat_cv$gamma_supp, stat_cv$beta_supp]
             option$gamma_supp = stat_cv$gamma_supp
           }else{
             result_cv = cv_all(X_1, y_1, D, option)
             beta_hat = result_cv$beta_hat
             stat_cv = result_cv$stat_cv
             X_new = X_2
             D_new = D
           }
           option$m = m
           option$beta_hat = beta_hat

           if (!is.null(option$nu)){
             nu_s = option$nu
           }else{
             nu_s = 10.^seq(0, 2, by = 0.2)
           }
           num_nu = length(nu_s)
           results = list()
           stats = list()

           for (i in 1:num_nu){
             nu = nu_s[i]
             stats[[i]] = sk.W_fixed(X_new, D_new, y_2, nu, option)
             W = stats[[i]]$W
             results[[i]] = list()
             method = 'knockoff'
             results[[i]]$sk = select(W, q, method)
             method = 'knockoff+'
             results[[i]]$sk_plus = select(W, q, method)
             stats[[i]]$nu = stat_cv$nu
             if (n_1 < p || (n-n_1) < (p+m)){
               stats[[i]]$gamma_supp = stat_cv$gamma_supp
             }
           }
         },
         'fill' = {
           if (!is.null(option$seed)){
             set.seed(option$seed)
           }
           rand_rank = sample(n)

           # record random order for function W_path
           option$rand_rank = rand_rank
           ind1 = rand_rank[1: floor(n * option$frac)]
           ind2 = rand_rank[floor(n*option$frac+1):length(rand_rank)]
           X_1 = X[ind1, ]
           y_1 = y[ind1]
           X_2 = X[ind2, ]
           y_2 = y[ind2]
           n_2 = nrow(X_2)

           if (p <= n_2 && n_2 < (p+m)){
             beta_hat = solve(t(X) %*% X) %*% t(X) %*% y
             sigma_hat = sqrt(sum((y - X %*% beta_hat)^2) / (n - p))
             add_number = p+m-n_2
             add_noise = stats::rnorm(add_number) * sigma_hat
             y_2 = c(y_2, add_noise)
             X_2 = rbind(X_2, matrix(0, add_number, p))
           }
           result_cv = cv_all(X_1, y_1, D, option)
           beta_hat = result_cv$beta_hat
           stat_cv = result_cv$stat_cv
           X_new = X_2
           D_new = D
           option$beta_hat = beta_hat
           option$m = m
           if (!is.null(option$nu)){
             nu_s = option$nu
           }else{
             nu_s = 10.^seq(0, 2, by = 0.2)
           }
           num_nu = length(nu_s)
           results = list()
           stats = list()
           for (i in 1:num_nu){
             nu = nu_s[i]
             stats[[i]] = sk.W_fixed(X_new, D_new, y_2, nu, option)
             W = stats[[i]]$W
             results[[i]] = list()
             method = 'knockoff'
             results[[i]]$sk = select(W, q, method)
             method = 'knockoff+'
             results[[i]]$sk_plus = select(W, q, method)
             stats[[i]]$nu = stat_cv$nu
           }
         },
         'appoint' = {
           option$beta_hat = option$beta_appointed
           if (!is.null(option$nu)){
             nu_s = option$nu
           }else{
             nu_s = 10.^seq(0, 2, by = 0.2)
           }
           num_nu = length(nu_s)
           results = list()
           stats = list()

           for (i in 1:num_nu){
             nu = nu_s[i]
             stats[[i]] = sk.W_fixed(X, D, y, nu, option)
             W = stats[[i]]$W
             results[[i]] = list()
             method = 'knockoff'
             results[[i]]$sk = select(W, q, method)
             method = 'knockoff+'
             results[[i]]$sk_plus = select(W, q, method)
           }
         })

  return(list(results = results, stats = stats))

}
