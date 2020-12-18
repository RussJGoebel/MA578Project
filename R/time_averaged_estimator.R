
#' Unbiased Time averaged estimator
#'
#' @param res A matrix returned by coupled_ARWM().
#' @param m
#' @param k
#' @param tau Tau as returned by coupled_ARWM(). Note that this is the first INDEX that the chain meets, so the meeting time is tau-1.
#' @param h A vectorized function.
#'
#' @return
#' @export
#'
#' @examples
time_averaged_estimator <- function(res,m,k,tau,h = function(x){x}){

  # CHECK that tau gets us the index at which they meet.
  # h should be vectorized.

  Res <- res$Res
  couple <- res$Res_couple

  biased_estimate <- 0
  for(l in ((k+1):(m+1))){
    biased_estimate <- biased_estimate + h(Res[l,])
  }

  biased_estimate <- biased_estimate/(m-k+1)

  bias_correction <- 0
  if((tau-1) >= k+2){
  for(l in (k+2):(tau-1)){
    bias_correction = bias_correction+min(1,(l-k)/(m-k+1))*(h(Res[l,])-h(couple[l-1,]))
  }
  }

  return(list(estimate = biased_estimate+bias_correction, biased_estimate = biased_estimate, bias_correction = bias_correction))

}
