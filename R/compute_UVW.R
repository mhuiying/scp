#' Internal function: compute U, V, W in Appendix A
#'
#' @param Q inverse of the covariance matrix
#' @param Y a vector with \eqn{n} data values.
#'
#' @return a list of U, V, W values
#' @export
compute_UVW <- function(Q,Y){
  Ysub     <- Y[-1]
  Qrow     <- Q[1,-1]
  Qsub     <- Q[-1,-1]
  dQsub    <- diag(Qsub)
  Yhat     <- -sum(Qrow*Ysub)/Q[1,1]
  Ytilde   <- -(Qsub%*%Ysub-dQsub*Ysub)/dQsub
  R        <- Ysub-Ytilde
  U        <- dQsub*R^2-Q[1,1]*Yhat^2
  V        <- Qrow*R + Q[1,1]*Yhat
  W        <- Qrow*Qrow/dQsub-Q[1,1]
  out      <- list(U=c(0,U),V=c(0,2*V),W=c(0,W))
  return(out)}
