#' Sample Maximal Coupling
#'
#' Uses maximal coupling to produce a single sample from distriubtions p and q.
#'
#'
#' @param psample A function which produces 1 sample from the distribution p
#' @param qsample A function 1 sample from the distribution q
#' @param pdensity A function of double x returning the log density of x from distribution p
#' @param qdensity A function of double x returning the log density of x from distribution q
#'
#' @return
#' @export
#'
#' @examples
#'
#' psamp <- function(){rnorm(1,0,1)}
#' qsamp <- function(){rexp(1,1)}
#' pden <- function(x){dnorm(x,0,1, log = TRUE)}
#' qden <- function(x){dexp(x,1, log = TRUE)}
#'
#'
#' niter <- 100000
#' rmat <- matrix(nrow = niter,ncol = 2)
#'
#' for(i in 1:niter){
#'
#' r <- sample_maximal_coupling(psamp,qsamp,pden,qden)
#' rmat[i,1] <- r$x
#' rmat[i,2] <- r$y
#'
#' }
#'
#' hist(rmat[,1])
#' hist(rmat[,2])
#'


sample_maximal_coupling <- function(psample,qsample,pdensity,qdensity){
  # pdensity and qdensity should be log densities

  x <- psample()
  y <- c()
  w <- log(runif(1,0,exp(pdensity(x))))

  if(w <= qdensity(x)){
    y <- x
  }
  else{

    run_condition <- T
    while(run_condition){
      yprop <- qsample()
      wstar <- log(runif(1,0,exp(qdensity(yprop))))
      if(wstar > pdensity(yprop)){
        y <- yprop
        run_condition = F
      }
    }


  }

  return(list(x = x, y = y))

}

#' Helper functions for maximal sampling
#'
#' Maximal sampling requires evaluating densities and sampling from distributions.
#' Given a multivariate normal distribution, these helper functions generate functions
#' to achieve this.
#'
#' @param mu Mean of distribution.
#' @param sigma Covariance matrix \eqn{\Sigma} to sample from
#' @param R Cholesky factor of \eqn{\Sigma}
#' @param log_scaling Paramater \eqn{\sigma} in \eqn{e^{2\sigma}\Sigma}
#'
#' @return
#' @export
#'
#' @examples
#'
#' m1 <- 10
#' sigma1 <- diag(1)
#' log_scaling1 <- 1
#'
#' m2 <- -10
#' sigma2 <- diag(1)
#' log_scaling2 <- 1
#'
#' h1 <- generate_mvn_sampler_functions(mu = m1,sigma  =sigma1, log_scaling = log_scaling1)
#' h2 <- generate_mvn_sampler_functions(mu = m2,sigma  =sigma2, log_scaling = log_scaling2)
#'
#' niter <- 10000
#' rmat <- matrix(nrow = niter,ncol = 2)
#'
#' for(i in 1:niter){
#'
#' r <- sample_maximal_coupling(h1$msample,h2$msample,h1$log_density,h2$log_density)
#' rmat[i,1] <- r$x
#' rmat[i,2] <- r$y
#'
#' }
#'
#' hist(c(rmat[,1],rmat[,2]))
#'
generate_mvn_sampler_functions <- function(mu = c(0,0),sigma = diag(length(mu)),R = t(chol(sigma)),log_scaling = 0){

  msample <- function(){sample_mvn(mu,sigma,R,log_scaling)}
  log_density <- function(x){log_density_mvn(x,mu,sigma,R,log_scaling)}

  return(list(msample = msample, log_density = log_density))

}





