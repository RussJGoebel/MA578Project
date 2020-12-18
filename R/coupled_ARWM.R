#' Coupled Adaptive Random Walk Metropolis Algorithm
#'
#' Given a log posterior density, generates samples using coupled adaptive random walk metropolis.
#' Code is modified from course materials from the Boston University course MA578 as taught by Yves Atchade.
#'
#' @param Niter Number of iterations for the algorithm to run
#' @param Y
#' @param X
#'
#' @return
#' @export
#' @examples
#'
#' x <- matrix(rnorm(100))
#' y <- rbinom(100,1,exp(x)/(1+exp(x)))
#'
#' logpi <- logpi_logistic
#'
#' r <- coupled_ARWM(1000,y,x)
#'
#'
#'
#'
coupled_ARWM = function(Niter,Y,X){
    p = length(X[1,]);
    alpha = 0.7;
    Res =matrix(NA,ncol=p,nrow=Niter);
    beta = numeric(p);
    lsig = 0;
    mu = beta;
    C = diag(p); R = t(chol(C))
    cpt_Update = 0; L_update = 500;
    lpi = logpi(beta,Y,X);

    chains_have_not_met <- T
    iterations_remain <- T
    jj <- 0
     while(iterations_remain){
      jj <- jj+1


      beta_prop = beta + exp(lsig)*as.vector(R%*%rnorm(p));
      lpi_prop = logpi(beta_prop,Y,X);


      A = exp(lpi_prop-lpi);
      Acc = min(1,A);
      if (runif(1)<=Acc){
        beta = beta_prop;
        lpi = lpi_prop;
      }
      Res[jj,] = beta;
      lsig = lsig + (1/jj^alpha)*(Acc -0.3);
      mu = mu + (1/jj)*(beta - mu);
      bmu = as.vector(beta-mu);
      C = C + (1/jj)*(outer(bmu,bmu) - C);
      if (cpt_Update == L_update){
        R = t(chol(C));
        cpt_Update = 0;
      }else{
        cpt_Update = cpt_Update + 1;
      }

      iterations_remain <- (jj < Niter)

     }

    return(Res)
}
