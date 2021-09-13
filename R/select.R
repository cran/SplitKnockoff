# R function for calculating the statistics W
# 1 June, 2021
# revised 7 July, 2021
# author：Haoxue Wang (haoxwang@student.ethz.ch)

#' splitknockoff.select
#'
#' @param W statistics W_j for testing null hypothesis
#' @param q target FDR
#'
#' @return S array of selected variable indices
#' @usage splitknockoff.select(W, q)
#' @export
#'
splitknockoff.select<-function(W, q){

  #  Inputs:
  #     W - statistics W_j for testing null hypothesis beta_j = 0.
  #       q - target FDR
  #
  #   Outputs:
  #       S - array of selected variable indices
  W <- t(W)
  t = sort(c(0,abs(W[W!=0])))
  ratio = matrix(0, 1, length(t))
  for (i in 1:length(t)){
  ratio[i] = sum(W <= -t[i]) / max(1, sum(W >= t[i]))
  }
  nindex <- which(ratio <= q & ratio!=0)
  if(is.null(nindex)==TRUE){
    T <- Inf
  }else{
  index <- nindex[1]
  T = t[index]
  }
  S = which(W >= T)
  structure(list(call = match.call(),
                 S = S),
            class = 'selectS_result')
}
