#' hitting point calculator on a given path
#' @description
#' calculate the hitting time and the sign of
#' respective variable in a path.
#'
#' @param coef the path for one variable
#' @param lambdas respective value of lambda in the path
#'
#' @return Z: the hitting time
#' @return r: the sign of respective variable at the hitting time
#' @export
#'
hittingpoint <- function(coef, lambdas){
  # hittingpoint calculates the hitting time and the sign of
  # respective variable in a path.
  #
  # input arguments
  # coef: the regularization path for one variable
  # lambdas: the respective values of lambda in the path
  #
  # output arguments
  # Z: the hitting time
  # r: the sign of respective variable at the hitting time

  n_lambda = length(lambdas)

  Z = 0
  r = 0

  # calculate Z and r
  for (j in 1:n_lambda){
    if (abs(coef[j]) != 0){
      Z = lambdas[j]
      r = sign(coef[j])
      break
    }
  }

  return(list(Z = Z, r = r))
}
