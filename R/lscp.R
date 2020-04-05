#' Local Spatial Conformal Prediction (LSCP) At a Single Input Location
#'
#' @description This function provides the local conformal prediction interval for spatial location \code{s0} given obserations \code{s, Y}.
#'
#' @param s0 prediction location
#' @param s an \eqn{n \times d}{n x d} matrix or data-frame with \eqn{d} coordinates of the \eqn{n} data locations.
#' @param Y a vector with \eqn{n} data values.
#' @param eta numerical value of the kernel bandwidth for the weight schema in conformal prediction. Defauls to \eqn{Inf} meaning equal weight on surrounding \eqn{m} points.
#' @param m an postive integer representing the number of nearest locations used for prediction. Depends on eta.
#' @param pred_fun
#' @param theta
#' @param alpha significance level. Defaults to 0.05.
#' @param dfun non-conformity measure
#' @param precision Defaults to 0.01.
#'
#' @return
#' @export
#'
#' @examples
lscp = function(s0,s,Y,eta=Inf,m=NULL,pred_fun=krige_pred,theta=NULL,alpha=0.05,
               dfun = "std_residual2",precision=0.01){

  if(identical(pred_fun, krige_pred)){
    if(dfun == "std_residual2"){

    } else if(dfun == "residual2"){

    } else{
      bins     <- seq(0.01,0.2,0.01)
      thetaHat <- get_theta(s,Y,dists=bins)
    }
  }

  idx = which(s[,1]==s0[1] & s[,2]==s0[2])
  if(length(idx) > 0){s = s[-idx,]; Y = Y[-idx]}

  dist = sqrt( (s0[1]-s[,1])^2 + (s0[2]-s[,2])^2 )
  these = (dist < 2*eta)
  dist_mat = apply(s[these,], 1, FUN = function(s0) sqrt( (s0[1]-s[,1])^2 + (s0[2]-s[,2])^2 ) )
  these = apply(dist_mat, 2, FUN = function(x) order(x)[1:15])
  these = unique(as.vector(these))
  if(is.null(m)) m = length(these)

  s     = rbind(s[these,], s0)
  Y     = c(Y[these], NA)
  d     = as.matrix(dist(s))
  if(eta==Inf){
    w <- rep(1/m,m)
  }else{
    w <- d[m+1,-(m+1)]/eta
    w <- exp(-0.5*w^2)
    w <- w/sum(w)
  }

  if(dfun == "std_residual2"){
    dfun = function(y,yhat,sd) (y-yhat)^2/sd^2
    need_sd = TRUE
  }else if(dfun == "residual2"){
    dfun = function(y,yhat,sd=1) (y-yhat)^2/sd^2
    need_sd = FALSE
  }else if(dfun == "abs_residual"){
    dfun = function(y,yhat,sd=1) abs(y-yhat)/sd
    need_sd = FALSE
  }else if(dfun == "abs_std_residual"){
    dfun = function(y,yhat,sd) abs(y-yhat)/sd
    need_sd = TRUE
  }else if(!is.function(dfun)){
    stop("Please provide a valid non-conformity measure.")
  }

  Y_cand = seq(min(Y[-(m+1)]),max(Y[-(m+1)]),precision)
  p_y = c()

  for(y in Y_cand){
    Y[m+1] = y
    delta = rep(NA, m+1)
    for(i in 1:(m+1)){
      if(need_sd)
        yi = pred_fun(as.numeric(s[i,]),s[-i,],Y[-i],alpha=alpha,thetaHat=thetaHat,return_sd=need_sd)
      else
        yi = list(yhat = pred_fun(as.numeric(s[i,]),s[-i,],Y[-i],alpha=alpha,thetaHat=thetaHat),
                  sd = 1)
      delta[i] = dfun(Y[i], yi$yhat, yi$sd)
    }
    p_y = c(p_y, sum(w*(delta >= delta[m+1])))
  }

  gamma = range(Y_cand[which(p_y >= floor((m+1) * alpha) / (m+1))])
  return(gamma)}
