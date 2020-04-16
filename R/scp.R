#' Spatial Conformal Prediction (SCP) At a Single Input Location
#'
#' @description This function provides the conformal prediction interval for spatial location \code{s0} given obserations \code{s, Y}.
#'
#' @param s0 prediction location
#' @param s an \eqn{n \times d}{n x d} matrix or data-frame with \eqn{d} coordinates of the \eqn{n} data locations.
#' @param Y a vector with \eqn{n} data values.
#' @param global logical; if TRUE , \code{scp} function returns the result of global spatial conformal prediction \code{\link{gscp}};
#' if FALSE, \code{scp} function returns the result of local spatial conformal prediction \code{\link{lscp}}.
#' @param eta numerical value of the kernel bandwidth for the weight schema in conformal prediction. Defauls to \eqn{Inf} meaning equal weight on surrounding \eqn{m} points.
#' @param m an postive integer representing the number of nearest locations used for prediction. Depends on eta.
#' @param pred_fun spatial point prediction function
#' @param alpha significance level. Defaults to 0.05.
#' @param dfun non-conformity measure
#' @param precision Defaults to 0.01.
#'
#' @return A vector of lower and upper bounds of the conformal prediction interval.
#' @export
#'
#' @examples
#' N = 41; n = N^2
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
#' idx = which(s[,1]==s0[1] & s[,2]==s0[2])
#' pred_fun = function(s0,s,Y,alpha=0.05) return(mean(Y))
#' PI = scp(s0,s[-idx,],Y[-idx],pred_fun=pred_fun, dfun="abs_residual")
#' cat(paste("True value: ", Y[idx], "\n"))
#' cat(paste("Prediction Interval: [ ", PI[1], ",", PI[2], "]"))

scp = function(s0,s,Y,global=FALSE,eta=Inf,m=NULL,pred_fun=krige_pred,alpha=0.05,
               dfun = "std_residual2",precision=NULL){

  idx = which(s[,1]==s0[1] & s[,2]==s0[2])
  if(length(idx) > 0){s = s[-idx,]; Y = Y[-idx]}

  if(!global){
    dist = sqrt( (s0[1]-s[,1])^2 + (s0[2]-s[,2])^2 )
    these = (dist < 2*eta)
    dist_mat = apply(s[these,], 1, FUN = function(s0) sqrt( (s0[1]-s[,1])^2 + (s0[2]-s[,2])^2 ) )
    these = apply(dist_mat, 2, FUN = function(x) order(x)[1:15])
    these = unique(as.vector(these))
    if(is.null(m)) m = length(these)
  } else{
    m = length(Y)
    these = 1:m
  }

  p_y = plausibity(s0, s[these,], Y[these], eta, m)

  gamma = range(Y_cand[which(p_y >= floor((m+1) * alpha) / (m+1))])
  return(gamma)

}
