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
#'
#'
#' r <- coupled_ARWM(1000,y,x, type = "reflection maximal", acceptance_probability = 0.3)
#'
#' n <- 1000
#' t <- numeric(n)
#' for(i in 1:n){
#' r <- coupled_ARWM(1,y,x, type = "reflection maximal")
#' t[i] <- r$tau
#'
#' if(i %% 10 == 0){print(i)}
#'
#' }
#'
#'
#'
coupled_ARWM = function(Niter,y,X,logpi = logpi_logistic, type = c("maximal","reflection maximal","none"),
                        acceptance_probability = 0.3,
                        beta_init = numeric(length(X[1,])),
                        beta_couple_init = numeric(length(X[1,])),
                        C_init= diag(length(X[1,]))){

    type <- match.arg(type)

    if(identical(type,"maximal")){

      coupling_sample <- function(ps,qs,R,log_scaling = lsig){

        beta_fns <- get_mvn_sampler_functions(mu = ps,R = R, log_scaling = lsig)
        beta_couple_fns <- get_mvn_sampler_functions(mu = qs,R = R,log_scaling = lsig)

        sample_maximal_coupling(beta_fns$msample,beta_couple_fns$msample,
                                beta_fns$log_density,beta_couple_fns$log_density)


      }

    }
    if(identical(type,"reflection maximal")){

      coupling_sample <- function(beta,beta_couple,R,log_scaling = lsig){
        sample_reflection_maximal(beta,beta_couple,sigma = tcrossprod(R)*exp(2*log_scaling))
      }

    }

    p = length(X[1,]);
    alpha = 0.7;
    Res = matrix(NA,ncol=p,nrow=Niter);
    Res_couple = Res


    beta <-  beta_init;
    beta_couple <- beta_couple_init;

    lpi = logpi(beta,y,X);
    lpi_couple = logpi(beta_couple,y,X);

    lsig = 1;
    mu = beta;
    C = C_init; R = t(chol(C))
    cpt_Update = 0; L_update = 500;



    chains_have_not_met <- T*!identical(type,"none")

    iterations_remain <- T
    jj <- 0
    tau <- 1


     while(iterations_remain | chains_have_not_met & jj < 20000){
      jj <- jj+1


      if(chains_have_not_met & (jj != 1)){ #(if jj == 1, we generate x1 before coupling).

       tau = tau+1

       # sample from maximal coupling ------------------------------------------------------------------
       maximal_sample_couple <- coupling_sample(beta,beta_couple,R,lsig)
       beta_prop <- maximal_sample_couple$x
       beta_couple_prop <- maximal_sample_couple$y

       # evaluate log densities
       lpi_prop <- logpi(beta_prop,y,X)
       lpi_couple_prop <- logpi(beta_couple_prop,y,X)

       # Obtain acceptance criterion
       A <- lpi_prop-lpi
       A_couple <- lpi_couple_prop-lpi_couple

       Acc <- min(0,A)
       Acc_couple <- min(0,A_couple)

       # Check acceptance against common uniform RV
       lu <- log(runif(1))
       if(lu <= Acc){
         beta = beta_prop;
         lpi = lpi_prop;
       }
       if(lu <= Acc_couple){
         beta_couple <-  beta_couple_prop;
         lpi_couple <- lpi_couple_prop
       }
       Res[jj,] = beta;
       Res_couple[jj,] = beta_couple

       if(isTRUE(all.equal(beta,beta_couple))){chains_have_not_met <- F}

      } else{

      (beta_prop = sample_mvn(mu = beta,R = R,log_scaling = lsig));
      lpi_prop = logpi(beta_prop,y,X);


      A = lpi_prop-lpi;
      Acc = min(0,A); # logged
      if (log(runif(1))<=Acc){
        (beta = beta_prop);
        lpi = lpi_prop;
      }
      Res[jj,] = beta;
      Res_couple[jj-1+(jj==1),] = ifelse(jj==1,NA,beta)
      }

      # update adaptive parameters -----------------------

      lsig = lsig + (1/jj^alpha)*(exp(Acc) - 0.3);
      mu = mu + (1/jj)*(beta - mu);
      bmu = as.vector(beta-mu);
      C = C + (1/jj)*(outer(bmu,bmu) - C)+C*(identical(bmu,numeric(p)))/jj;
      if (cpt_Update == L_update){
        R = t(chol(C));
        cpt_Update = 0;
      }else{
        cpt_Update = cpt_Update + 1;
      }

      iterations_remain <- (jj < Niter)

      if((jj == dim(Res)[1]) & (chains_have_not_met)){

        rownum <- dim(Res)[1]
        colnum <- dim(Res)[2]

        Res <- rbind(Res,matrix(NA,nrow = rownum,ncol = colnum))
        Res_couple <- rbind(Res_couple,matrix(NA,nrow = rownum,ncol = colnum))


      }
      chains_have_not_met #DELETE DEBUG

     }

    remove_tail_na <- complete.cases(Res)

    if(identical(type,"none")) Res_couple <- NULL

    return(list(Res = Res[remove_tail_na,], Res_couple = Res_couple[remove_tail_na,],tau = tau, C = C))
}
