#' split knockoff selector given W statistics
#'
#' @param W statistics W_j for testing null hypothesis
#' @param q target FDR
#' @param method option$method can be 'knockoff' or 'knockoff+'
#'
#' @return S: array of selected variable indices
#' @export
#'
select <- function(W, q, method = 'knockoff+'){
  #   Inputs:
  #       W - statistics W_j for testing null hypothesis beta_j = 0.
  #       q - target FDR
  #       method - either 'knockoff' or 'knockoff+'
  #                Default: 'knockoff+'
  #
  #   Outputs:
  #       S - array of selected variable indices

  T = threshold(W, q, method)
  S = which(W >= T)
  return(S)
}
