#' compute the threshold for variable selection
#'
#' @param W statistics W_j for testing null hypothesis beta_j = 0
#' @param q target FDR
#' @param method option$method can be 'knockoff' or 'knockoff+'
#'
#' @return T: threshold for variable selection
#' @export
#'
threshold <- function(W, q, method = 'knockoff+'){
  #   Inputs:
  #       W - statistics W_j for testing null hypothesis beta_j = 0.
  #       q - target FDR
  #       method - either 'knockoff' or 'knockoff+'
  #                Default: 'knockoff'
  #
  #   Outputs:
  #       T - threshold for variable selection

  if (method == 'knockoff')
    offset = 0
  else if (method == 'knockoff+')
    offset = 1
  else
    stop("Invalid threshold method")

  t = sort(c(0,abs(W[W!=0])))
  ratio = rep(0, length(t))
  for (i in 1:length(t)){
    ratio[i] = (offset + sum(W <= -t[i])) / max(1, sum(W >= t[i]))
  }
  nindex <- which(ratio <= q)
  if (is.null(nindex) == TRUE){
    T <- Inf
  }else{
    index <- nindex[1]
    T = t[index]
  }

  return(T)
}
