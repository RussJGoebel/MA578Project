#' Log Posterior Density for Logistic Regression
#'
#' Evaluates the log posterior density for logistic regression using a particular prior.
#' Code is modified from course materials from the Boston University course MA578 as taught by Yves Atchade.
#'
#' @param beta Beta values
#' @param y Binary vector on which to fit logistic regression
#' @param X Matrix of continous predictors on which to fit logistic regression
#' @param c_prior Paramater for prior.
#'
#' @return
#' @export
#'
#' @examples
#'
#' x <- matrix(rnorm(100))
#' y <- rbinom(100,1,exp(x)/(1+exp(x)))
#' logpi_logistic(0,y,x)
#'
logpi_logistic = function(beta,y,X,prior = 100){
  #log-density for logistic regression

  mp = X%*%beta;
  val = sum(y*mp -log(1+exp(mp)))
  val = val -(0.5/prior^2)*crossprod(beta);
  return(val)
}


#' Log Posterior Density for Linear Regression
#'
#' Evaluates the log posterior density for linear regression using a particular prior.
#' Code is modified from course materials from the Boston University course MA578 as taught by Yves Atchade.
#' @param beta
#' @param prior
#'
#' @return
#' @export
#'
#' @examples
#'
#' x <- matrix(rnorm(100))
#' x <- cbind(1,x)
#' y <- x%*%c(5,1)+rnorm(50,0,1/2)
#' logpi_lm(c(5,1),y,x)
#'
#'
logpi_lm <- function(beta = c(0,0),y,X, prior = 10^(-3)*numeric(length(beta))){

  -(1/2)*(sum((y-X%*%beta)^2)+crossprod(prior*beta,beta))

}
