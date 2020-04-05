#' Kriging Prediction Function
#'
#' @param s0 prediction location
#' @param s an \eqn{n \times d}{n x d} matrix or data-frame with \eqn{d} coordinates of the \eqn{n} data locations.
#' @param Y a vector with \eqn{n} data values.
#' @param alpha significance level. Defaults to 0.05.
#' @param thetaHat spatial covariance parameters as in \code{\link{mat_cov}}.
#' If not given, emperical variogram is used to estimate \code{thetaHat}.
#' @param interval logical; if TRUE, \code{Krige_pred} function returns prediction interval;
#' if FALSE, \code{Krige_pred} function returns point prediction interval. Defaults to FALSE.
#' @param return_sd logical; if TRUE, \code{Krige_pred} function returns
#' standard deviation along with the point prediction. Defaults to FALSE.
#'
#' @return a value of point prediction if \code{interval} is FALSE
#' or a vector of lower and upper bounds of the Kriging prediction intervalif \eqn{interval} is TRUE.
#' @export
#'
#' @examples
#' N = 21; n = N^2
#' S = seq(0,1,length=N)
#' s = expand.grid(S,S)
#' d = as.matrix(dist(s))
#'
#' theta        = c(0,3,0.1,0.7)
#' names(theta) = c("Nugget","PartialSill","Range","Smoothness")
#' C = mat_cov(d,theta)
#' X = t(chol(C))%*%rnorm(n)
#' Y = X^3 + rnorm(n)
#'
#' s0 = c(0.5, 0.5)
#' krige_pred(s0,s,Y)
#' krige_pred(s0,s,Y,interval=TRUE)
krige_pred = function(s0,s,Y,alpha=0.05,thetaHat=NULL,interval=FALSE,return_sd=FALSE){

  # Suppose s does not have duplicates
  # check s0 is numeric
  # check s is data.frame. What if s is in matrix format

  if( length(Y) != nrow(s) )
    stop( paste("the number of Y obs,", length(Y), ", does not match the number of locations,", nrow(s)) )
  if( length(s0) != ncol(s) )
    stop( paste("the dimension of s0,", length(s0), ", does not match the dimension of s,", ncol(s)) )

  # Estimate Matern covariance parameters if not given
  if(is.null(thetaHat)){
    bins     <- seq(0.01,0.2,0.01)
    thetaHat <- get_theta(s,Y,dists=bins)
  }

  idx = which(s[,1]==s0[1] & s[,2]==s0[2])
  if(length(idx) > 0){s = s[-idx,]; Y = Y[-idx]}

  # distance matrix
  s_all = rbind(s0,s)
  d_all = as.matrix(dist(s_all))

  # inverse covariance matrix
  Q = solve(mat_cov(d_all,thetaHat))

  # Kriging prediction
  yhat  = as.numeric(-Q[1,-1]%*%Y/Q[1,1])
  sd = 1/sqrt(Q[1,1])
  gamma = yhat + qnorm(1-alpha/2)*c(-1,1)*sd

  yhat = ifelse(return_sd, list(yhat=yhat,sd=sd), yhat)
  ifelse(interval, return(gamma), return(yhat) )

}
