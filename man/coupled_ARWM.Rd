% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coupled_ARWM.R
\name{coupled_ARWM}
\alias{coupled_ARWM}
\title{Coupled Adaptive Random Walk Metropolis Algorithm}
\usage{
coupled_ARWM(Niter, Y, X)
}
\arguments{
\item{Niter}{Number of iterations for the algorithm to run}

\item{X}{}
}
\value{

}
\description{
Given a log posterior density, generates samples using coupled adaptive random walk metropolis.
Code is modified from course materials from the Boston University course MA578 as taught by Yves Atchade.
}
\examples{

x <- matrix(rnorm(100))
y <- rbinom(100,1,exp(x)/(1+exp(x)))

logpi <- logpi_logistic

r <- coupled_ARWM(1000,y,x)




}