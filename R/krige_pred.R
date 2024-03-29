#' Kriging prediction function
#'
#' @description This function provides an example for \code{pred_fun} in \code{\link{scp}},
#' \code{\link{plausibility}}, and \code{\link{plausibility_contour}}, which provides a point prediction
#' for location \code{s0} (and corresponding standard error), given obserations \code{s} and \code{Y}.
#'
#' @param s0 prediction location, a numeric vector with \code{length = 2}.
#' @param s an \eqn{n \times 2}{n x 2} \code{matrix} or \code{data.frame} with two coordinates of \eqn{n} locations.
#' @param Y a vector with \eqn{n} values corresponding to \code{Y(s)}.
#' @param return_sd logical. if \code{TRUE}, \code{Krige_pred} function returns the
#' standard error of \code{Y(s0)} along with the point prediction;
#' if \code{FALSE}, \code{Krige_pred} function only returns the point prediction. Defaults to \code{FALSE}.
#'
#' @return If \code{return_sd} is FALSE (default), the output is a value of point prediction for \code{Y(s0)};
#'         If \code{return_sd} is TRUE, the output is a \code{list} with the following elements:
#'    \item{yhat }{point prediction for \code{Y(s0)}}
#'    \item{sd }{standard error for \code{Y(s0)}}
#'
#' @export
#'
#' @examples
#' #?sample_data
#' s0 = c(0.5,0.5)
#' s  = sample_data$s
#' Y  = sample_data$Y
#'
#' krige_pred(s0,s,Y)
#' krige_pred(s0,s,Y,return_sd=TRUE)
krige_pred = function(s0,s,Y,return_sd=FALSE){

  if( length(Y) != nrow(s) )
    stop( paste("the number of Y obs,", length(Y), ", does not match the number of locations,", nrow(s)) )
  if( length(s0) != ncol(s) )
    stop( paste("the dimension of s0,", length(s0), ", does not match the dimension of s,", ncol(s)) )

  bins     = seq(0.01,0.2,0.01)
  thetaHat = get_theta(s,Y,dists=bins)

  idx = which(s[,1]==s0[1] & s[,2]==s0[2])
  if(length(idx) > 0){s = s[-idx,]; Y = Y[-idx]}

  # distance matrix
  s_aug = rbind(s0,s)
  d_aug = as.matrix(dist(s_aug))

  # inverse covariance matrix
  Q = solve(mat_cov(d_aug,thetaHat))

  # Kriging prediction
  yhat  = as.numeric(-Q[1,-1]%*%Y/Q[1,1])
  sd    = 1/sqrt(Q[1,1])

  if(!return_sd)
    return(yhat)
  else
    return(list(yhat=yhat,sd=sd))

}
