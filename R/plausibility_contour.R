#' generate plausibility contour
#'
#' @description This function provides the plausibility contour for \code{Y(s0)},
#' given observations \code{s} and \code{Y},
#' using spatial conformal prediction algorithms.
#'
#' @param s0 prediction location, a numeric vector with \code{length = 2}.
#' @param s an \eqn{n \times 2}{n x 2} \code{matrix} or \code{data.frame} with two coordinates of \eqn{n} locations.
#' @param Y a vector with \eqn{n} values corresponding to \code{Y(s)}.
#' @param global logical; if \code{TRUE} , \code{scp} function returns the result of global spatial conformal prediction (GSCP);
#' if \code{FALSE}, \code{scp} function returns the result of local spatial conformal prediction (LSCP) and users need to
#' specify \code{eta < Inf} or \code{m} \eqn{\le} \code{n}. Defaults to \code{TRUE}.
#' @param eta kernel bandwidth for weight schema, a positve scalar with smaller value meaning more localized procedure.
#' Defauls to \code{Inf}, which puts equal weight on surrounding \code{m} points.
#' @param m an postive integer representing the number of nearest locations to use for prediction.
#' Default to \code{NULL}. If \code{global = TRUE}, \code{m = n};
#' if \code{global = FALSE} and \code{m} is not specified, \code{m} would be determined by \code{eta}.
#' @param pred_fun spatial prediction function with inputs being \code{s0, s, Y} and ouputs being predicted \code{Y(s0)}
#' (and its standard error). Defaults to \code{\link{krige_pred}}.
#' @param thetaHat a vector of Matern parameters, representing nugget, partial sill, range, and smoothness as in Mao. et al. (2020).
#'             Defaults to \code{NULL}. It will be ignored if \code{pred_fun} is not \code{krige_pred}.
#' @param dfun non-conformity measure with four options.
#'             In which, \code{"residual2"} (default) represents squared residual and
#'             \code{"std_residual2"} represents standardized squared residual.
#' @param precision a positive scalar represents how dense \code{Y(s)} candidates (\code{Y_cand}) are. Defaults to \code{NULL}.
#'
#' @return The output is a \code{data.frame} of \code{Y_cand} and corresponding plausibility values.
#' @export
#'
#' @author Huiying Mao, \email{hmao@@samsi.info}, Brian Reich \email{bjreich@@ncsu.edu}
#' @references to be entered
#' @seealso \code{\link{plausibility}}
#'
#' @examples
#' ## generate plausibility contour for Y(s0), where s0 = c(0.5,0.5), using sample data
#'
#' #?sample_data
#' s0 = c(0.5,0.5)
#' s  = sample_data$s
#' Y  = sample_data$Y
#'
#' pc = plausibility_contour(s0=s0,s=s,Y=Y)
#' plot(pc)
#'
#' idx = which(s[, 1] == s0[1] & s[, 2] == s0[2])
#' abline(v = Y[idx], col = "red", lty = 2)
#' legend("topright", col=1:2, lty=1:2, c("plausibility", "true value"))
#'
plausibility_contour = function(s0,s,Y,global=TRUE,eta=Inf,m=NULL,pred_fun=krige_pred,thetaHat=NULL,
                                dfun=c("residual2","std_residual2"),precision=NULL){
  dfun = match.arg(dfun)
  .prime(s0,s,Y,global,eta,m,dfun)
  Y_cand = .generate_Y_cand(pred_fun, dfun, precision)
  p_df = .plausibility_contour(Y_cand,s_aug,Y_aug,w_aug,d_aug,pred_fun,thetaHat,T_dfun)
  p_df = new_plausibility_contour(p_df)
  return(p_df)

}

#' internal plausibility_contour calculation function
#'
#' @param Y_cand a vector of candidate Y values.
#' @param s_aug augmented locations, as an output of \code{\link{.prime}}.
#' @param Y_aug augmented Y values, as an output of \code{\link{.prime}}.
#' @param w_aug augmented weights, as an output of \code{\link{.prime}}.
#' @param pred_fun spatial prediction function with inputs being \code{s0, s, Y} and ouputs being predicted \code{Y(s0)}
#' (and its standard error). Defaults to \code{\link{krige_pred}}.
#' @param T_dfun transformed \code{dfun}, as an output of \code{\link{trans_dfun}}.
#'
#' @return a data frame of Y_cand and p_y.
#' @export
#' @keywords internal

.plausibility_contour = function(Y_cand,s_aug,Y_aug,w_aug,d_aug,pred_fun,thetaHat,T_dfun){

  dfun = T_dfun$fun

  if(identical(pred_fun, krige_pred)){
    if(is.null(thetaHat)){
        bins     = seq(0.01,0.2,0.01)
        thetaHat = get_theta(s,Y,dists=bins)
    }
    Q = solve(mat_cov(d_aug,thetaHat))

    if( T_dfun$residual2 )
      p_df = .fast_plausibility_contour(Y_aug,w_aug,Q,std=T_dfun$std)
    else
      p_df = .Q_plausibility_contour(Q,Y_cand,Y_aug,w_aug,pred_fun,T_dfun)
  }else{
    n = length(Y_aug)
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
#' @keywords internal
#'
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
#' @param Y_aug augmented Y values, as an output of \code{\link{.prime}}.
#' @param Y_cand a vector of candidate Y values.
#' @param Q an \eqn{n+1 \times n+1} inverse of estimated covariance matrix.
#' @param Y_aug augmented Y values, as an output of \code{\link{.prime}}.
#' @param w_aug augmented weights, as an output of \code{\link{.prime}}.
#' @param pred_fun spatial prediction function with inputs being \code{s0, s, Y} and ouputs being predicted \code{Y(s0)}
#' (and its standard error). Defaults to \code{\link{krige_pred}}.
#' @param T_dfun transformed \code{dfun}, as an output of \code{\link{trans_dfun}}.
#'
#' @return
#' @export
#' @keywords internal
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

new_plausibility_contour = function(x = list()){
  stopifnot(is.list(x))
  structure(x, class = "plausibility_contour")
}

#' Plot a plausibility_contour object
#'
#' @param df a plausibility_contour object, returned by \code{\link{plausibility_contour}}
#' @param ...
#'
#' @export
#' @keywords internal
plot.plausibility_contour = function(df, ...){
  plot(df[[1]], df[[2]], type = "l",
       xlab = "Y candidates", ylab = "plausibility",
       ...)
}
