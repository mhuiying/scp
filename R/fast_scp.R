#' Fast Spatial Conformal Prediction Intervals
#'
#' @description Internal function provides the conformal prediction interval
#' for spatial location \code{s0} given partial obserations \code{s1,...,sm, Y1,.., Ym}
#' when the square of the standard kriging residuals are used as the non-conformality measures,
#' and weights are provided.
#'
#' @param s0 prediction location
#' @param s an \eqn{n \times d}{n x d} matrix or data-frame with \eqn{d} coordinates of the \eqn{n} data locations.
#' @param Y a vector with \eqn{n} data values.
#' @param thetaHat estimated Matern covariance parameters as in \code{\link{mat_cov}}
#' @param w weights for the non-conformity measures.
#' @param alpha significance level. Defaults to 0.05.
#'
#' @return A vector of lower and upper bounds of the conformal prediction interval.
#' @export

fast_scp = function(s0,s,Y,thetaHat,alpha=0.05){

  if( length(Y) != nrow(s) )
    stop( paste("the number of Y obs,", length(Y), ", does not match the number of locations,", nrow(s)) )
  if( length(s0) != ncol(s) )
    stop( paste("the dimension of s0,", length(s0), ", does not match the dimension of s,", ncol(s)) )
  if( length(w) != nrow(s) )
    stop( paste("the number of weights,", length(w), ", does not match the number of locations,",nrow(s)) )
  if( sum(w)!=1 )
    stop(paste("Please double check weights. The summation is,", sum(w), ", NOT 1"))

  m = length(Y)
  s = rbind(s,s0)
  Y = c(Y,NA)
  d = as.matrix(dist(s))

  Q   = solve(mat_cov(d,thetaHat))
  UVW = compute_UVW(Q,Y)
  U   = UVW$U[-(m+1)]
  V   = UVW$V[-(m+1)]
  W   = UVW$W[-(m+1)]

  low    = (-V + sqrt(V^2 - 4*U*W))/(2*W)
  upp    = (-V - sqrt(V^2 - 4*U*W))/(2*W)
  lu_w   = c(w,w)
  lu_lab = c(rep(1, length(low)), rep(-1, length(upp)))
  lu_w   = lu_w[order(c(low, upp))]
  lu_lab = lu_lab[order(c(low, upp))]

  y0    = sort(c(low, upp))
  p_y   = cumsum(lu_w*lu_lab)
  gamma = range(y0[which(p_y >= floor((m+1) * alpha) / (m+1))])
  return(gamma)
}
