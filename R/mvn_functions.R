#' Multivariate Normal Sample
#'
#' Function to generate a single sample from the multivariate normal distribution.
#' Allows easy parameterization of the covariance matrix as \eqn{e^{2\sigma}\Sigma}.
#'
#' @param mu Mean of distribution to sample from
#' @param sigma Covariance matrix \eqn{\Sigma} to sample from
#' @param R Cholesky factor of \eqn{\Sigma}
#' @param log_scaling Paramater \eqn{\sigma} in \eqn{e^{2\sigma}\Sigma}
#'
#' @return
#' @export
#'
#' @example
#'
#' n <- 1000
#' s1 <- numeric(n)
#' s2 <- numeric(n)
#'
#' mu <- c(1,1)
#' sigma <- diag(2)
#' log_scaling = 1
#'
#' for(i in 1:n){
#'  r <- sample_mvn(mu,sigma,log_scaling = log_scaling)
#'  s1[i] <- r[1]
#'  s2[i] <- r[2]
#' }
#'
#' plot(s1~s2)
#'
#'
sample_mvn <- function(mu,
                       sigma = diag(length(mu)),
                       R = t(chol(sigma)),
                       log_scaling = 0){

  p <- length(mu)
  mu + exp(log_scaling)*as.vector(R%*%rnorm(p))

}


#' Multivariate Normal Log Density
#'
#' Function to evaluate the logged multivariate normal density at a point.
#'
#' @param x Point to evaluate density at.
#' @param mu Mean of distribution.
#' @param sigma Covariance matrix \eqn{\Sigma} to sample from
#' @param R Cholesky factor of \eqn{\Sigma}
#' @param log_scaling Paramater \eqn{\sigma} in \eqn{e^{2\sigma}\Sigma}
#'
#' @return
#' @export
#'
#'
log_density_mvn <- function(x,
                            mu = numeric(length(x)),
                            sigma = diag(length(mu)),
                            R = t(chol(sigma)),
                            log_scaling = 0){

  p <- length(x)

  # Compute Quadratic Form
  a <- x-mu
  q <- sum(forwardsolve(R,a)^2)*exp(-2*log_scaling)

  # Compute Determinant of Sigma
  d <- prod(diag(R))

  result <- -1/2*(q+p*log(2*pi)+p*2*log_scaling)-log(d)
  return(result)

}
