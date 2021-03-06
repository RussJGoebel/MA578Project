% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvn_functions.R
\name{log_density_mvn}
\alias{log_density_mvn}
\title{Multivariate Normal Log Density}
\usage{
log_density_mvn(
  x,
  mu = numeric(length(x)),
  sigma = diag(length(mu)),
  R = t(chol(sigma)),
  log_scaling = 0
)
}
\arguments{
\item{x}{Point to evaluate density at.}

\item{mu}{Mean of distribution.}

\item{sigma}{Covariance matrix \eqn{\Sigma} to sample from}

\item{R}{Cholesky factor of \eqn{\Sigma}}

\item{log_scaling}{Paramater \eqn{\sigma} in \eqn{e^{2\sigma}\Sigma}}
}
\value{

}
\description{
Function to evaluate the logged multivariate normal density at a point.
}
