#' plausibility_contour calculation
#'
#' @param Y_cand
#' @param s0
#' @param s
#' @param Y
#' @param global
#' @param eta
#' @param m
#' @param pred_fun
#' @param dfun
#'
#' @return
#' @export
#'
#' @examples
plausibility_contour = function(Y_cand=NULL,s0,s,Y,global,eta,m,pred_fun,dfun){

  T_dfun = trans_dfun(dfun)

  if(is.null(Y_cand))
    if( !identical(pred_fun, krige_pred) | !T_dfun$residual2 )
      Y_cand = seq(min(Y),max(Y),length.out = 100)

  sYw_aug = .hyper(s0,s,Y,global,eta,m)
  return(.plausibility_contour(Y_cand,sYw_aug,pred_fun,T_dfun))

}

#' internal plausibility_contour calculation function
#'
#' @param Y_cand
#' @param sYw_aug
#' @param pred_fun
#' @param T_dfun
#'
#' @return a data frame of Y_cand and p_y
#' @export
#' @keywords internal
#'
#' @examples
.plausibility_contour = function(Y_cand,sYw_aug,pred_fun,T_dfun){

  s_aug = sYw_aug$s_aug
  Y_aug = sYw_aug$Y_aug
  w_aug = sYw_aug$w_aug
  d_aug = as.matrix(dist(s_aug))
  dfun = T_dfun$fun

  if(identical(pred_fun, krige_pred)){
    bins     = seq(0.01,0.2,0.01)
    thetaHat = get_theta(s,Y,dists=bins)
    Q        = solve(mat_cov(d_aug,thetaHat))

    if( T_dfun$residual2 )
      p_df = .fast_plausibility_contour(Y_aug,w_aug,Q,std=T_dfun$std)
    else
      p_df = .Q_plausibility_contour(Q,Y_cand,Y_aug,w_aug,pred_fun,T_dfun)
  }else{
    n    = length(Y_aug)
    p_y = c()
    for(y in Y_cand){
      Y_aug[1] = y
      Yhat_aug = rep(NA,n)
      sd_aug   = rep(1 ,n)
      for(i in 1:n)
        Yhat_aug[i] = pred_fun(as.numeric(s_aug[i,]),s_aug[-i,],Y_aug[-i])
      if(T_dfun$std)
        delta = dfun(Y_aug, Yhat_aug, sd=sd_aug)
      else
        delta = dfun(Y_aug, Yhat_aug, sd=1)
      p_y = c(p_y, sum(w_aug*(delta >= delta[1])))
    }
    p_df = data.frame(Y_cand=Y_cand,p_y=p_y)
  }

  return(p_df)
}

#' Fast plausibility_contour calculation when pred_fun = Krige_fun and (dfun = "residual2" or "std_residual2")
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

.fast_plausibility_contour = function(Y_aug,w_aug,Q,std=TRUE){

  if(std)
    UVW = .compute_UVW(Q,Y_aug)
  else
    UVW = .compute_UVW2(Q,Y_aug)
  U   = UVW$U[-1]
  V   = UVW$V[-1]
  W   = UVW$W[-1]

  low    = (-V + sqrt(V^2 - 4*U*W))/(2*W)
  upp    = (-V - sqrt(V^2 - 4*U*W))/(2*W)

  w      = w_aug[-1]
  lu_w   = c(w,w)
  lu_lab = c(rep(1, length(low)), rep(-1, length(upp)))
  lu_w   = lu_w[order(c(low, upp))]
  lu_lab = lu_lab[order(c(low, upp))]

  y0    = sort(c(low, upp))
  p_y   = cumsum(lu_w*lu_lab)

  return(data.frame(Y_cand=y0,p_y=p_y))
}

#' Need inverse covariance matrix to calculate plausibility_contour
#'
#' @param Y_aug
#' @param w_aug
#' @param Q
#' @param std
#'
#' @return
#' @export
#'
#' @examples
.Q_plausibility_contour = function(Q,Y_cand,Y_aug,w_aug,pred_fun,T_dfun){

  dfun = T_dfun$fun
  p_y = c()
  for(y in Y_cand){
    Y_aug[1] = y
    Yhat_aug = as.numeric(-Q%*%Y_aug/diag(Q))
    sd_aug   = as.numeric(sqrt(diag(Q)))
    if(T_dfun$std)
      delta = dfun(Y_aug, Yhat_aug, sd=sd_aug)
    else
      delta = dfun(Y_aug, Yhat_aug, sd=1)
    p_y = c(p_y, sum(w_aug*(delta >= delta[1])))
  }
  return(data.frame(Y_cand=Y_cand,p_y=p_y))

}
