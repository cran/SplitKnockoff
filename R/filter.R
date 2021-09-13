# R function for filter of structural sparsity problem
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)

#' splitknockoff.filter
#'
#' Split Knockoff filter for structural sparsity problem.
#'
#' @param X the design matrix
#' @param D the linear transform
#' @param y the response vector
#' @param option options for creating the Knockoff statistics
#' option$eta: the choice of eta for creating the knockoff copy
#' option$q: the desired FDR control bound
#' option$method: 'knockoff' or 'knockoff+'
#' option$stage0: choose the method to conduct split knockoff.
#'       'fixed': fixed intercept assignment for PATH ORDER method.
#'          option$beta : the choice of fixed beta for step 0:
#'               'mle': maximum likelihood estimator.
#'               'ridge': ridge regression choice beta with lambda = 1/nu.
#'               'cv_split': cross validation choice of split LASSO over nu
#'               and lambda.
#'               'cv_ridge': cross validation choice of ridge regression
#'               over lambda.
#'       'path': take the regularization path of split LASSO as the intercept
#'               assignment for PATH ORDER method.
#'       'magnitude': using MAGNITUDE method.
#' option$lambda: a set of lambda appointed for path calculation
#' option$nu: a set of nu used for Split Knockoffs
#' option$normalize: whether to normalize the data
#'
#' @return results: a cell with the selected variable set in each cell w.r.t. nu.
#' @return Z: a cell with the feature significance Z in each cell w.r.t. nu.
#' @return t_Z: a cell with the knockoff significance tilde_Z in each cell w.r.t. nu.
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
#' filter_result <- splitknockoff.filter(X, D, y, option)
#' Z_path <- filter_result$Z
#' t_Z_path <- filter_result$t_Z
#' @export
#'
splitknockoff.filter <- function(X, D, y, option)
# Input Argument
# X : the design matrix.
# y : the response vector.
# D : the linear transform.
# option: options for creating the Knockoff statistics.
#	option$eta : the choice of eta for creating the knockoff copy.
#	option$q: the desired FDR control bound.
#	option$method: 'knockoff' or 'knockoff+'.
#	option$stage0: choose the method to conduct split knockoff.
#       'fixed': fixed intercept assignment for PATH ORDER method.
#           option.beta : the choice of fixed beta for step 0:
#               'mle': maximum likelihood estimator.
#               'ridge': ridge regression choice beta with lambda = 1/nu.
#               'cv_split': cross validation choice of split LASSO over nu
#               and lambda.
#               'cv_ridge': cross validation choice of ridge regression
#               over lambda.
#       'path': take the regularization path of split LASSO as the intercept
#               assignment for PATH ORDER method.
#       'magnitude': using MAGNITUDE method.
#	option$lambda: a set of lambda appointed for path calculation.
#	option$nu: a set of nu used for Split Knockoffs.
#	option$normalize: whether to normalize the data.

# Output Argument
# results: a cell with the selected variable set in each cell w.r.t. nu.
# Z: a cell with the feature significance Z in each cell w.r.t. nu.
# t_Z: a cell with the knockoff significance tilde_Z in each cell w.r.t. nu.
{
  if(option$normalize == "true"){
    X <- normc(X) # normalize(X)
    y <- normc(y) # normalize(y)
  }
  nu_s = option$nu
  num_nu <- length(nu_s)
  n_nu <- length(nu_s)
  Z <- list(data=NA,dim=c(num_nu, 1,n_nu))
  t_Z <- list(data=NA,dim=c(num_nu, 1,n_nu))

  n<- nrow(X)
  q = option$q
  method = option$method
  stage0 = option$stage0


  results <- list(data=NA,dim=c(num_nu, 1,n_nu))
  # if(all.equal(option$stage0, 'fixed') && all.equal(option$beta, 'ridge'))
  # {
  #   split_knockoffs.statistics.pathorder.fixed_beta(X, y, D, option)
  #   option$beta_choice <- beta_choice
  #
  # }
  for (i in 1: n_nu) {
    nu = nu_s[i]
    # filter.choose(option$stage0)
    path_result <- split_knockoffs.statistics.pathorder.W_path(X, D, y, nu, option)
    W <- path_result$W
    Z[[i]] <- path_result$Z
    t_Z[[i]] <- path_result$t_Z
    select_results <- splitknockoff.select(W, q)
    result <- select_results$S
    results[[i]]<-result
  }

  structure(list(call = match.call(),
                 results = results,
                 W = W,
                 Z = Z,
                 t_Z = t_Z),
            class = 'filter_result')
}


#' cv_filter
#'
#' Split Knockoff filter for structural sparsity problem, using cross validation
#'
#' @param X the design matrix
#' @param D the response vector
#' @param y the linear transform
#' @param option options for creating the Knockoff statistics
#' option$eta  the choice of eta for creating the knockoff copy
#' option$q the desired FDR control bound
#' option$method 'knockoff' or 'knockoff+'
#' option$stage0 choose the method to conduct split knockoff
#'       'fixed': fixed intercept assignment for PATH ORDER method.
#'          option$beta : the choice of fixed beta for step 0:
#'               'mle': maximum likelihood estimator.
#'               'ridge': ridge regression choice beta with lambda = 1/nu.
#'               'cv_split': cross validation choice of split LASSO over nu
#'               and lambda.
#'               'cv_ridge': cross validation choice of ridge regression
#'               over lambda.
#'       'path': take the regularization path of split LASSO as the intercept
#'               assignment for PATH ORDER method.
#'       'magnitude': using MAGNITUDE method.
#' option$lambda: a set of lambda appointed for path calculation
#' option$nu: a set of nu used for Split Knockoffs
#' option$k_fold: the fold used in cross validation
#' option$cv_rule: the rule used in CV
#'       'min': choose nu with minimal CV loss.
#'       'complexity': choose nu with minimal model complexity in the range
#'        of 0.99 * CV_loss <= min(CV_loss).
#'
#' @return result: selected features of Split Knockoffs with CV optimal selection of nu.
#' @return CV_loss: the CV loss of Split Knockoffs w.r.t. nu.
#' @return nu_optimal: the CV optimal nu.
#' @usage cv_filter(X, D, y, option)
#' @export
#'
#'
cv_filter <- function(X, D, y, option)
  #  Input Argument
  # X : the design matrix.
  # y : the response vector.
  # D : the linear transform.
  # option: options for creating the Knockoff statistics.
  #	option$eta : the choice of eta for creating the knockoff copy.
  #	option$q: the desired FDR control bound.
  #	option$method: 'knockoff' or 'knockoff+'.
  #	option$stage0: choose the method to conduct split knockoff.
  #       'fixed': fixed intercept assignment for PATH ORDER method.
  #           option$beta : the choice of fixed beta for step 0:
#               'mle': maximum likelihood estimator.
#               'ridge': ridge regression choice beta with lambda = 1/nu.
#               'cv_split': cross validation choice of split LASSO over nu
#               and lambda.
#               'cv_ridge': cross validation choice of ridge regression
#               over lambda.
#       'path': take the regularization path of split LASSO as the intercept
#           assignment for PATH ORDER method.
#       'magnitude': using MAGNITUDE method.
#	option$lambda: a set of lambda appointed for path calculation.
#	option$nu: a set of nu used for Split Knockoffs.
# 	option$k_fold: the fold used in cross validation.
# 	option$cv_rule: the rule used in CV.
#       'min': choose nu with minimal CV loss.
#       'complexity': choose nu with minimal model complexity in the range
#           of 0.99 * CV_loss <= min(CV_loss).
#
# Output Argument
# result: selected features of Split Knockoffs with CV optimal selection of nu.
# CV_loss: the CV loss of Split Knockoffs w.r.t. nu.
# nu_optimal: the CV optimal nu.
{
  nu <- option$nu
  nu_s <- option$nu
  loss_result<- split_knockoffs.statistics.pathorder.W_path(X, D, y, nu, option)
  results <-loss_result$results
  CV_loss <- loss_result$CV_loss
  complexities <- loss_result$complexities

  switch(option$cv_rule,
         min = {
           index <- which(CV_loss == min(CV_loss))
           nu_optimal <- nu_s(index)
           result <- results[index]},

         complexity= {
           set = which(0.99 * CV_loss <= min(CV_loss))
           nu_potential = nu_s[set]
           complexities = complexities[set]
           result_cv = results[set]
           index = which(complexities == min(complexities))
           nu_optimal = nu_potential[index]
           result = result_cv[index]
         }
  )

  structure(list(call = match.call(),
                 nu_optimal = nu_optimal,
                 CV_loss = CV_loss,
                 result = result),
            class = 'cv_filter_result')
}





