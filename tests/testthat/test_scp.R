context("testing the scp functionality")

test_that("Whether scp gives the right output", {

  source("tests/testthat/setup.R")

  # situation 1
  # kriging + residual2
  PI1 = scp(s0=s0,s=s[-idx,],Y=Y[-idx])
  expect_equal(PI1, c(-22.2214369049998, 10.992719070945))

  # situation 2
  # kriging + std_residual2
  PI2 = scp(s0=s0,s=s[-idx,],Y=Y[-idx],global=FALSE,eta=Inf,m=NULL,pred_fun=krige_pred,dfun="std_residual2")
  expect_equal(PI2, c(-22.2214369049998, 10.992719070945))

  # situation 3
  # kriging(Q) + abs_residual
  PI3 = scp(s0=s0,s=s[-idx,],Y=Y[-idx],global=FALSE,eta=Inf,m=NULL,pred_fun=krige_pred,dfun="abs_residual")
  expect_equal(PI3, c(-22.5475033874506, 20.4126287790842))

})
