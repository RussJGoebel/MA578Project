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


#' Sample using Reflection Maximal Sampling
#'
#' See Unbiased Markov chain Monte Carlo with Couplings, by Pierre E. Jacob, John O'leary, and Yves F. Atchade for details.
#'
#'
#' @param Sigma Covariance matrix of distributions p and q.
#' @param mu_p Mean of the distribution p.
#' @param mu_q Mean of the distribution q.
#' @param s_sample A function which returns one sample from the distribution s. Defaults to standard multivariate normal.
#' @param s_density A function of x returning the log density of x from the distribution s. Defaults to standard multivariate normal.
#'
#' @return
#' @export
#'
#' @examples
#'
#' mu_p = -c(1,1)
#' mu_q = c(1,1)
#' sigma = diag(2)
#'
#'
#'
#'
#'
#' niter <- 1000
#' rmat <- matrix(nrow = niter, ncol = 4)
#'
#' for(i in 1:niter){
#' r <- sample_reflection_maximal(mu_p,mu_q,sigma)
#' rmat[i,1:2] <- r$x
#' rmat[i,3:4] <- r$y
#' }
#'
#' hist(c(rmat[,1],rmat[,2]))
#' plot(rmat[,3]~rmat[,4])
#'
#'
#'
sample_reflection_maximal <- function(mu_p,mu_q,sigma,
                                      s_sample = function(){sample_mvn(numeric(length(mu_p)))},
                                      s_density = function(x){log_density_mvn(x)}){



    if(length(mu_p) > 1){

      qdq <- eigen(sigma, symmetric = TRUE)
      q_factor<- qdq$vectors

      root_sigma_inverse <- q_factor%*% diag(1/sqrt(qdq$values)) %*% t(q_factor)
      root_sigma <- q_factor%*% diag(sqrt(qdq$values)) %*% t(q_factor)
    } else{
      root_sigma_inverse <- 1/sqrt(sigma)
      root_sigma <- sqrt(sigma)
    }


    x <- s_sample()

    z <- root_sigma_inverse %*% (mu_p-mu_q)
    e <- z/sqrt(sum(z^2))

    u <- runif(1)

    A <- min(1,exp(s_density(x+z)-s_density(x)))

    y <- x+z
    if(u > A){
      y <- x-c(2*(crossprod(e,x))) * e
    }

    x <- mu_p+root_sigma %*% x
    y <- mu_q+root_sigma %*% y

    x <- round(x,digits = 6)
    y <- round(y,digits = 6)



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
#' h1 <- get_mvn_sampler_functions(mu = m1,sigma  =sigma1, log_scaling = log_scaling1)
#' h2 <- get_mvn_sampler_functions(mu = m2,sigma  =sigma2, log_scaling = log_scaling2)
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
get_mvn_sampler_functions <- function(mu = c(0,0),sigma = diag(length(mu)),R = t(chol(sigma)),log_scaling = 0){

  msample <- function(){sample_mvn(mu,sigma,R,log_scaling)}
  log_density <- function(x){log_density_mvn(x,mu,sigma,R,log_scaling)}

  return(list(msample = msample, log_density = log_density))

}



