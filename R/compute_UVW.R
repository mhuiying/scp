#' compute U, V, W
#'
#' @description compute U, V, W in Mao et al. (2020) Appendix A when delta_i = (Y_i - Yhat_i)^2 / v_i^2
#'
#' @param Q inverse of the covariance matrix
#' @param Y a vector with \eqn{n} data values.
#'
#' @return a list of U, V, W values
#' @keywords internal
#'
#' @export
#'
.compute_UVW = function(Q,Y){
  Ysub     = Y[-1]
  Qrow     = Q[1,-1]
  Qsub     = Q[-1,-1]
  dQsub    = diag(Qsub)
  Yhat     = -sum(Qrow*Ysub)/Q[1,1]
  Ytilde   = -(Qsub%*%Ysub-dQsub*Ysub)/dQsub
  R        = Ysub-Ytilde
  U        = dQsub*R^2-Q[1,1]*Yhat^2
  V        = Qrow*R + Q[1,1]*Yhat
  W        = Qrow*Qrow/dQsub-Q[1,1]
  out      = list(U=c(0,U),V=c(0,2*V),W=c(0,W))
  return(out)}

#' compute U, V, W
#' @description compute U, V, W in Mao et al. (2020) Appendix A when delta_i = (Y_i - Yhat_i)^2
#'
#' @param Q inverse of the covariance matrix
#' @param Y a vector with \eqn{n} data values.
#'
#' @return a list of U, V, W values
#' @keywords internal
#'
.compute_UVW2 = function(Q,Y){
  Ysub     = Y[-1]
  Qrow     = Q[1,-1]
  Qsub     = Q[-1,-1]
  dQsub    = diag(Qsub)
  Yhat     = sum(Qrow*Ysub)/Q[1,1]
  Ytilde   = Qsub%*%Ysub/dQsub
  U        = Ytilde^2-Yhat^2
  V        = Ytilde*Qrow/dQsub - Yhat
  W        = Qrow*Qrow/(dQsub*dQsub)-1
  out      = list(U=c(0,U),V=c(0,2*V),W=c(0,W))
  return(out)}
