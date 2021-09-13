# R function for simulation
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)

# this function calculate the FDR and Power when gamma^* is available, i.e.
# in simulations.

#' simu_eval
#'
#' calculate the FDR and Power when gamma^* is available, i.e. in simulations.
#'
#' @param gamma_true true signal of gamma
#' @param result the estimated support set of gamma
#'
#' @return fdr: false discovery rate of the estimated support set
#' @return power: power of the estimated support set
#' @usage
#' simu_eval(gamma_true, result)
#' @export
#'
simu_eval <- function(gamma_true, result)
# input argument
# gamma_true: true signal of gamma
# result: the estimated support set of gamma

# output argument
# fdr: false discovery rate of the estimated support set
# power: power of the estimated support set
{
    total_posi <- length(result)
    false_posi <- total_posi - sum(gamma_true[result] != 0)
    fdr = false_posi / total_posi
    power <- sum(gamma_true[result] != 0) / sum(gamma_true != 0)
  structure(list(call = match.call(),
                 fdr = fdr,
                 power = power),
            class = 'eval_result')

}



#' simu_unit_cv
#'
#' the simulation unit for simulation experiments with cross validation.
#'
#' @param n the sample size
#' @param p the dimension of variables
#' @param D_s the set of linear transform
#' @param A SNR
#' @param c feature correlation
#' @param k number of nonnulls in beta
#' @param option option for split knockoffs
#'
#' @return simu_data: a structure contains the following elements
#' @return simu_data$fdr: fdr of cv optimal nu in split knockoffs
#' @return simu_data$power: power of cv optimal nu in split knockoffs
#' @return simu_data$fdr_knock: fdr of knockoffs
#' @return simu_data$cv_list: length(D_s) * length(nu_s) matrix with cv loss for
#' each nu for each D
#' @return simu_data$chosen_nu: cv selected nu for each D
#' @import MASS
#' @usage
#' simu_unit_cv (n, p, D_s, A, c, k, option)
#' @export
#'
#'
simu_unit_cv <- function(n, p, D_s, A, c, k, option)
# input argument
# n: the sample size
# p: the dimension of variables
# D_s: the set of linear transform
# A: SNR
# c: feature correlation
# k: number of nonnulls in beta
# option: option for split knockoffs

# output argument
# simu_data: a structure contains the following elements
# simu_data$fdr: fdr of cv optimal nu in split knockoffs
# simu_data$power: power of cv optimal nu in split knockoffs
# simu_data$fdr_knock: fdr of knockoffs
# simu_data$power_knock: power of knockoffs
# simu_data$cv_list: length(D_s) * length(nu_s) matrix with cv loss for
# each nu for each D
# simu_data$chosen_nu: cv selected nu for each D

{
  sigma = 1 # noise level
  test <- option$test
  num_nu <- length(option$nu)

  # generate X
  Sigma = matrix(0, p, p)
  for (i in 1: p) {
    for (j in 1: p) {
      Sigma[i, j] = c^(abs(i - j))
    }
  }
  set.seed(test)
  X <- mvrnorm(matrix(0, p, 1), Sigma, n) # package MASS needed
  # generate beta and gamma
  beta_true <- matrix(0,p, 1)
  for (i in 1: k) {
    beta_true[i, 1] = A
    if(i == 1){
      beta_true[i, 1] = -A
    }
  }
  # create matrices to store results
  fdr_cv <- matrix(0, 3, 1)
  power_cv <- matrix(0, 3, 1)
  cv_list <- matrix(0, 3, num_nu)
  chosen_nu <- matrix(0, 3, 1)
  fdr_knockoff <- matrix(0, 2, 1)
  power_knockoff <- matrix(0, 2, 1)

  ############# begin simulation ############
  for (D_choice in 1: 3) {
    D = D_s[D_choice]

    m <- nrow(D)
    gamma_true <- D %*% beta_true

    # generate varepsilon
    set.seed(test)

    # generate noise and y
    varepsilon <- matrix(0,n, 1) %*% sqrt(sigma)
    y <- X %*% beta_true + varepsilon


    # running knockoff as a comparison
    # if (m <= p){
    # result <- split_knockoffs.private.convert_knockoff(X, D, y, option)
    # list[fdr_knockoff(D_choice, 1), power_knockoff(D_choice, 1)] <- split_knockoffs.private.simu_eval(gamma_true, result)
    # }

    cv_result <- cv_filter(X, D, y, option)
    result <- cv_result$result
    cv_list[D_choice, ] <- cv_result$cv_list
    chosen_nu[D_choice] <- cv_result$chosen_nu
    eval_result <- simu_eval(gamma_true, result)
    fdr_cv[D_choice, 1] <- eval_result$fdr
    power_cv[D_choice, 1] <- eval_result$power
  }
  # compute the means

  structure(list(call = match.call(),
                 fdr_cv = fdr_cv,
                 power_cv = power_cv,
                 cv_list =cv_list,
                 chosen_nu=chosen_nu),
            class = 'cv_data')


}



#' simu_unit
#'
#' the simulation unit for simulation experiments.
#'
#' @param n the sample size
#' @param p the dimension of variables
#' @param D the linear transform
#' @param A SNR
#' @param c feature correlation
#' @param k number of nonnulls in beta
#' @param option option for split knockoffs
#'
#' @return simu_data: a structure contains the following elements
#' @return simu_data$fdr_split: a vector recording fdr of split knockoffs w.r.t.nu
#' @return simu_data$power_split: a vector recording power of split knockoffs w.r.t.nu
#' @return simu_data$fdr_knock: fdr of knockoffs
#' @return simu_data$power_knock: power of knockoffs
#' @usage
#' simu_unit(n, p, D, A, c, k, option)
#' @export
#'
#'
simu_unit <- function(n, p, D, A, c, k, option)
# input argument
# n: the sample size
# p: the dimension of variables
# D: the linear transform
# A: SNR
# c: feature correlation
# k: number of nonnulls in beta
# option: option for split knockoffs

# output argument
# simu_data: a structure contains the following elements
# simu_data$fdr_split: a vector recording fdr of split knockoffs w.r.t.nu
# simu_data$power_split: a vector recording power of split knockoffs w.r.t.nu
# simu_data$fdr_knock: fdr of knockoffs
# simu_data$power_knock: power of knockoffs
{
  sigma = 1 # noise level
  tests = 20 # number of experiments
  num_nu <- length(option$nu)
  m = nrow(D)

  # generate X
  Sigma = matrix(0,p, p)
  for (i in 1: p) {
    for (j in 1: p) {
      Sigma[i, j] = c^(abs(i - j))
    }
  }
  set.seed(test)
  X <- mvrnorm(n, matrix(0, p, 1), Sigma) # package MASS needed
  # generate beta and gamma
  beta_true <- matrix(0,p,1)
  for (i in 1: k) {
    beta_true[i, 1] <- A
    if(i %% 3 == 1){
      beta_true[i, 1] <- -A
    }
  }
  gamma_true <- D %*% beta_true

  # create matrices to store results
  fdr_split <- matrix(0,tests, num_nu)
  power_split <- matrix(0,tests, num_nu)
  #fdr_knockoff <- matrix(0,tests, 1)
  #power_knockoff <- matrix(0,tests, 1)

  ########### begin simulation ##########
  for (test in 1: tests) {
    # generate varepsilon
    set.seed(test)

    # generate noise and y
    varepsilon <- rnorm(n) * sqrt(sigma)
    y <- X %*% beta_true + varepsilon

    # if(m <= p){
    #   result <- split_knockoffs.private.convert_knockoff(X, D, y, option)
    #   simu_eval(gamma_true, result)
    #   fdr_knockoff[test]<-fdr
    #   power_knockoff[test] <- power
    #    }
    filter_result <- splitknockoff.filter(X, D, y, option)
    results <- filter_result$results

    for (i in 1: num_nu) {
      result <- results[[i]]
      eval_result <- simu_eval(gamma_true, result)
      fdr_split[test,i]<-eval_result$fdr
      power_split[test,i] <- eval_result$power
    }
  }
  # compute the means

  mean_fdr_split <- apply(fdr_split,2,mean)
  mean_power_split <- apply(power_split,2,mean)

  # mean_fdr_knockoff <- mean(fdr_knockoff)
  # mean_power_knockoff <- mean(power_knockoff)

  structure(list(call = match.call(),
                 fdr_split = mean_fdr_split,
                 power_split = mean_power_split
                 # fdr_knock = fdr_knockoff,
                 # power_knock = power_knockoff
                 ),
            class = 'simu_data')

}
