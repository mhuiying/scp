#' Fast Spatial Conformal Prediction Intervals
#'
#' @description Internal function provides the conformal prediction interval
#' for spatial location \code{s0} given partial obserations \code{s1,...,sm, Y1,.., Ym}
#' when the square of the standard kriging residuals are used as the non-conformality measures,
#' and weights are provided.
#'
#' @param s0 prediction location
#' @param s an \eqn{n \times d} matrix or data-frame with \eqn{d} coordinates of the \eqn{n} data locations.
#' @param Y a vector with \eqn{n} data values.
#' @param w weights for the non-conformity measures.
#' @param Q an \eqn{n+1 \times n+1} inverse of estimated covariance matrix
#' @param alpha significance level. Defaults to 0.05.
#'
#' @return A vector of lower and upper bounds of the conformal prediction interval.
#' @export

fast_scp = function(s0,s,Y,w,Q,alpha=0.05,std=TRUE){

  if( length(Y) != nrow(s) )
    stop( paste("the number of Y obs,", length(Y), ", does not match the number of locations,", nrow(s)) )
  if( length(s0) != ncol(s) )
    stop( paste("the dimension of s0,", length(s0), ", does not match the dimension of s,", ncol(s)) )
  if( length(w) != nrow(s) )
    stop( paste("the number of weights,", length(w), ", does not match the number of locations,",nrow(s)) )
  if( sum(w)!=1 )
    stop(paste("Please double check weights. The summation is,", sum(w), ", NOT 1"))
  if( any(dim(Q) != nrow(s)+1) )
    stop( paste("the dimension of Q,", dim(Q)[1],"by", dim(Q)[2],
                ", does not match the number of augmented locations,",nrow(s)+1) )

  Y_aug = c(NA,Y)
  if(std)
    UVW = .compute_UVW(Q,Y_aug)
  else
    UVW = .compute_UVW2(Q,Y_aug)
  U   = UVW$U[-1]
  V   = UVW$V[-1]
  W   = UVW$W[-1]

  low    = (-V + sqrt(V^2 - 4*U*W))/(2*W)
  upp    = (-V - sqrt(V^2 - 4*U*W))/(2*W)
  lu_w   = c(w,w)
  lu_lab = c(rep(1, length(low)), rep(-1, length(upp)))
  lu_w   = lu_w[order(c(low, upp))]
  lu_lab = lu_lab[order(c(low, upp))]

  y0    = sort(c(low, upp))
  p_y   = cumsum(lu_w*lu_lab)

  m     = length(Y)
  gamma = range(y0[which(p_y >= floor((m+1) * alpha) / (m+1))])
  return(gamma)
}
