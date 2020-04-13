context("testing the conformal_pred functionality")

test_that("Whether fast_scp gives the right output", {

  set.seed(123)
  N = 11; n = N^2
  S = seq(0,1,length=N)
  s = expand.grid(S,S)
  d = as.matrix(dist(s))

  theta        = c(0,3,0.1,0.7)
  names(theta) = c("Nugget","PartialSill","Range","Smoothness")
  C = mat_cov(d,theta)
  X = t(chol(C))%*%rnorm(n)
  Y = X^3 + rnorm(n)

  # Estimate spatial covariance parameters
  bins     = seq(0.01,0.2,0.01)
  thetaHat = get_theta(s,Y,dists=bins)

  # spatial prediction
  s0  = c(0.5, 0.5)
  idx = which(s[,1]==s0[1] & s[,2]==s0[2])
  PI  = conformal_pred(s0,s[-idx,],Y[-idx],thetaHat,eta=0.1)

  expect_equal(PI, c(-21.061147745455, 16.6843427225051))
})
