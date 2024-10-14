#' calculate the CV optimal beta
#' @description
#' cv_all calculate the CV optimal beta
#' in the problem 1/n |y - X beta|^2 + 1/nu |D beta - gamma|^2 + lambda |gamma|_1.
#'
#' @param X the design matrix
#' @param y the response vector
#' @param D the linear transform
#' @param option options for screening
#'
#' @return beta_hat: CV optimal beta
#' @return stat_cv: various intermedia statistics
#' @export
#'
#' @importFrom glmnet glmnet
cv_all <- function(X, y, D, option){
  # cv_all calculate the CV optimal beta
  # in the problem 1/n |y - X beta|^2 + 1/nu |D beta - gamma|^2 + lambda |gamma|_1.
  #
  # input arguments
  # X : the design matrix
  # y : the response vector
  # D : the linear transform
  #
  # output arguments
  # beta_hat: CV optimal beta
  # stat_cv: various intermedia statistics

  k_fold = 5
  n = nrow(X)
  m = nrow(D)
  p = ncol(D)

  # appoint the set of \nu
  if (!is.null(option$nu_cv)){
    nu_s = option$nu_cv
  }else{
    nu_s = 10.^seq(0, 2, by = 0.4)
  }

  # appoint a set of lambda
  if (!is.null(option$lambda_cv)){
    lambda_s = option$lambda_cv
  }else{
    lambda_s = 10.^seq(0, -8, by = -0.4)
  }
  nlambda = length(lambda_s)

  test_size = floor(n / k_fold)

  # randomly split data

  set.seed(1)
  rand_rank = sample(n)

  # create matrix to store result for split lasso test loss
  loss_sl = array(0, dim = c(k_fold, length(nu_s), nlambda))

  for (k in 1:k_fold){
    # generate test set
    test_index = test_size * (k-1) + (1: test_size)
    test = rand_rank[test_index]

    # training set
    X_train = X[!(1:n %in% test), ]
    y_train = y[!(1:n %in% test)]

    # test set
    X_test = X[test, ]
    y_test = y[test]

    # normalization
    n_train = nrow(X_train)
    n_test = nrow(X_test)
    X_train = X_train - matrix(rep(colMeans(X_train), n_train), nrow = n_train, ncol = p, byrow = TRUE)
    y_train = y_train - mean(y_train)
    X_test = X_test - matrix(rep(colMeans(X_test), n_test), nrow = n_test, ncol = p, byrow = TRUE)
    y_test = y_test - mean(y_test)

    # test loss for split lasso
    for (i in 1:length(nu_s)){
      nu = nu_s[i]

      # generate split LASSO design matrix
      A_beta <- rbind(X_train/sqrt(n_train),D/sqrt(nu))
      A_gamma <- rbind(matrix(0,n_train,m),-diag(m)/sqrt(nu))
      tilde_y <- rbind(matrix(y_train, ncol = 1)/sqrt(n_train),matrix(0,m,1))

      # set penalty
      penalty <- matrix(1,m+p,1)

      for (tmp in 1: p) {
        penalty[tmp, 1] = 0
      }

      # lasso path settings for glmnet

      fit = glmnet(cbind(A_beta,A_gamma) , tilde_y, lambda = lambda_s, standardize = FALSE, intercept = FALSE, penalty.factor = penalty)
      coefs <- fit$beta
      for (j in 1:length(lambda_s)){
        coef = coefs[ ,j]
        beta = coef[1:p]
        # calculate loss
        y_sl = as.matrix(X_test) %*% as.matrix(beta)
        loss_sl[k, i, j] = norm(y_sl - y_test, type = "F")^2/test_size
      }
    }
  }

  mean_loss_sl = colMeans(loss_sl)

  # find minimal
  nu_number = which(mean_loss_sl == min(mean_loss_sl), arr.ind = TRUE)[1, 1]
  lambda_number = which(mean_loss_sl == min(mean_loss_sl), arr.ind = TRUE)[1, 2]
  nu_sl = nu_s[nu_number]
  lambda_sl = lambda_s[lambda_number]

  stat_cv = list()
  stat_cv$nu = nu_sl
  stat_cv$lambda = lambda_sl

  # calculate beta
  A_beta <- rbind(X/sqrt(n),D/sqrt(nu_sl))
  A_gamma <- rbind(matrix(0,n,m),-diag(m)/sqrt(nu_sl))
  tilde_y <- rbind(matrix(y, ncol = 1)/sqrt(n),matrix(0,m,1))

  fit = glmnet(cbind(A_beta,A_gamma) , tilde_y, lambda = lambda_sl, standardize = FALSE, intercept = FALSE, penalty.factor = penalty)
  coef <- fit$beta
  beta_hat = coef[1:p]

  return(list(beta_hat = beta_hat, stat_cv = stat_cv))
}
