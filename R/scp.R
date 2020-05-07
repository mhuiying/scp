#' Spatial conformal prediction at input location(s)
#'
#' @description This function provides the spatial conformal prediction interval for location(s) \code{s0},
#'  given obserations \code{s} and \code{Y}.
#'
#' @param s0 prediction location(s), a numeric vector with \code{length = 2},
#'                                   or a \code{matrix} with \code{ncol = 2},
#'                                   or a \code{data.frame} with two coordinates.
#' @param s an \eqn{n \times 2}{n x 2} \code{matrix} or a \code{data.frame} with two coordinates of \eqn{n} locations.
#' @param Y a vector with \eqn{n} values corresponding to \code{Y(s)}.
#' @param global logical; if \code{TRUE} , \code{scp} function returns the result of global spatial conformal prediction (GSCP);
#' if \code{FALSE}, \code{scp} function returns the result of local spatial conformal prediction (LSCP)
#' and users need to specify \code{eta}. Defaults to \code{TRUE}.
#' @param eta kernel bandwidth for weight schema, a positve scalar with smaller value meaning more localized procedure.
#' Defauls to \code{Inf}, which puts equal weight on surrounding \eqn{m} points.
#' @param m an postive integer representing the number of nearest locations to use for prediction.
#' Default depands on \code{eta}.
#' @param pred_fun spatial prediction function with inputs being \eqn{s0, s, Y} and ouputs being predicted \code{Y(s0)}
#' (and its standard error). Defaults to \code{\link{krige_pred}} representing Kriging prediction.
#' @param thetaHat a vector of Matern parameters, representing nugget, partial sill, range, and smoothness as in Mao. et al. (2020).
#'             Defaults to \code{NULL}. It will be ignored if \code{pred_fun} is not \code{krige_pred}.
#' @param dfun non-conformity measure with four options.
#'             In which, \code{"residual2"} (default) represents squared residual and
#'             \code{"std_residual2"} represents standardized squared residual.
#' @param precision a positive scalar represents how dense the candidates for \code{Y(s)} are. Defaults to \code{NULL}.
#' @param alpha significance level. Defaults to 0.05.
#'
#' @return The output is a \code{data.frame} of lower and upper bounds of the conformal prediction interval(s) for \code{s0}.
#' @export
#'
#' @author Huiying Mao, \email{hmao@@samsi.info}, Brian Reich \email{bjreich@@ncsu.edu}
#' @references to be entered
#' @seealso \code{\link{plausibility}}, \code{\link{plausibility_contour}}
#'
#' @examples
#' ## generate prediction interval for prediction locations(s) s0(s) using sample data
#'
#' #?sample_data
#' s  = sample_data$s
#' Y  = sample_data$Y
#'
#' # locations to predict
#' s0  = c(0.5,0.5)
#' s0s = rbind(c(0.4, 0.4), c(0.5,0.5), c(0.6, 0.6))
#'
#' # default prediction interval
#' scp(s0=s0,s=s,Y=Y)
#' scp(s0=s0s,s=s,Y=Y)
#'
#' # user define eta=0.1, where LSCP is considered
#' scp(s0=s0,s=s,Y=Y,eta=0.1)
#'
#' # user define non-conformity measure
#' scp(s0=s0,s=s,Y=Y,dfun="std_residual2")
#'
#' # user define prediction function
#' fun = function(s0,s,Y) return(mean(Y))
#' scp(s0=s0,s=s,Y=Y,pred_fun=fun)

scp = function(s0,s,Y,global=TRUE,eta=Inf,m=NULL,pred_fun=krige_pred,thetaHat=NULL,
               dfun=c("residual2","std_residual2"),precision=NULL,alpha=0.05){

  dfun = match.arg(dfun)

  if( !is.matrix(s0) & !is.data.frame(s0) ){
    s0 = matrix(s0, ncol = 2)
    colnames(s0) = colnames(s)
  }

  n_pred_loc = nrow(s0)
  gamma = data.frame(lower = rep(NA, n_pred_loc), upper = rep(NA, n_pred_loc))
  for(ith_loc in 1:n_pred_loc){
    .prime(s0[ith_loc,],s,Y,global,eta,m,dfun)
    Y_cand = .generate_Y_cand(pred_fun, dfun, precision)

    p_df   = .plausibility_contour(Y_cand,s_aug,Y_aug,w_aug,d_aug,pred_fun,thetaHat,T_dfun)
    Y_cand = p_df$Y_cand
    p_y    = p_df$p_y

    gamma[ith_loc,] = range(Y_cand[which(p_y >= floor((M+1) * alpha) / (M+1))])
  }

  return(gamma)

}
