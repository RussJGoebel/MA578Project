


test_that("log density function works", {
  library(mvtnorm)
  set.seed(1)
  lsig <- 3
  x <- c(1,2,3)
  Sigma <- matrix(rnorm(9),nrow = 3)
  Sigma <- crossprod(Sigma)
  package_value <- mvtnorm::dmvnorm(x,mean = c(1,2,3),sigma = exp(lsig)*Sigma, log = TRUE)
  myfun_value <- log_density_mvn(x,mu = c(1,2,3),sigma = Sigma,log_scaling = lsig)

  expect_equal(round(package_value,5),round(myfun_value,5))

})

