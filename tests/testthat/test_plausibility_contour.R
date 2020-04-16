context("testing the .plausibility_contour functionality")

test_that("Whether .plausibility_contour gives the right output", {

  source("tests/testthat/setup.R")

  # situation 1
  # kriging + residual2
  p_df1 = plausibility_contour(s0=s0,s=s[-idx,],Y=Y[-idx],global=FALSE,eta=Inf,m=NULL,pred_fun=krige_pred,dfun="residual2")
  expect_equal(p_df1[c(1,nrow(p_df1)),],
               structure(list(Y_cand = c(-58.9977770941897, 47.779809781105),
                              p_y = c(0.00826446280991736, 0)), row.names = c(1L, 240L), class = "data.frame"))

  # situation 2
  # kriging + std_residual2
  p_df2 = plausibility_contour(s0=s0,s=s[-idx,],Y=Y[-idx],global=FALSE,eta=Inf,m=NULL,pred_fun=krige_pred,dfun="std_residual2")
  expect_equal(p_df2[c(1,nrow(p_df2)),],
               structure(list(Y_cand = c(-58.929980325642, 47.712014664932),
                              p_y = c(0.00826446280991736, 0)), row.names = c(1L, 240L), class = "data.frame"))

  # situation 3
  # kriging(Q) + abs_residual
  p_df3 = plausibility_contour(s0=s0,s=s[-idx,],Y=Y[-idx],global=FALSE,eta=Inf,m=NULL,pred_fun=krige_pred,dfun="abs_residual")
  expect_equal(p_df3[c(1,nrow(p_df3)),],
               structure(list(Y_cand = c(-22.5475033874506, 60.8456943475876),
                              p_y = c(0.0578512396694215, 0.00826446280991736)),
                         row.names = c(1L, 100L), class = "data.frame"))
})
