#' calculate plausibility for \code{Y0}
#'
#' @description This function provides the plausibility of \code{Y(s0)} being \code{Y0},
#' given observations \code{s} and \code{Y},
#' using spatial conformal prediction algorithms.
#'
#' @param Y0 a scalar or a vector
#' @param s0 prediction location, a numeric vector with \code{length = 2}.
#' @param s an \eqn{n \times 2}{n x 2} \code{matrix} or \code{data.frame} with two coordinates of \eqn{n} locations.
#' @param Y a vector with \eqn{n} values corresponding to \code{Y(s)}.
#' @param global logical; if \code{TRUE} , \code{scp} function returns the result of global spatial conformal prediction (GSCP);
#' if \code{FALSE}, \code{scp} function returns the result of local spatial conformal prediction (LSCP) and users need to
#' specify \code{eta < Inf} or \code{m} \eqn{\leq} \code{n}. Defaults to \code{TRUE}.
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
#'             In which, \code{"residual2"} (default) represents squared residual,
#'             \code{"std_residual2"} represents standardized squared residual,
#'             \code{"abs_residual"} represents absolute residual,
#'             and \code{"std_abs_residual"} represents standardized absolute residual.
#'
#' @return The output is a scalar or a vector with plausibility values for \code{Y0}. The numbers are between 0 and 1.
#' @export
#'
#' @author Huiying Mao, \email{hmao@@samsi.info}, Brian Reich \email{bjreich@@ncsu.edu}
#' @references to be entered
#' @seealso \code{\link{plausibility_contour}}
#'
#' @examples
#' ## To predict Y(s0), where s0 = c(0.5,0.5), using sample data
#' ## What's the plausibility if Y(s0) = 0? Y(s0) = 1.5?
#'
#' #?sample_data
#' s0 = c(0.5,0.5)
#' s  = sample_data$s
#' Y  = sample_data$Y
#'
#' # plausibility for Y(s0) = 0: 0.8744795
#' plausibility(Y0=0,s0=s0,s=s,Y=Y)
#'
#' # plausibility for Y(s0) = 1.5: 0.4669839
#' plausibility(Y0=1.5,s0=s0,s=s,Y=Y)
#'
#' # plausibility for a sequence of Y0's
#' plausibility(Y0=seq(0,1,0.1),s0=s0,s=s,Y=Y)
#'
plausibility = function(Y0,s0,s,Y,global=TRUE,eta=Inf,m=NULL,pred_fun=krige_pred,thetaHat=NULL,
                        dfun=c("residual2","std_residual2")){

  dfun = match.arg(dfun)
  .prime(s0,s,Y,global,eta,m,dfun)
  if( !identical(pred_fun, krige_pred) | !grepl("residual2", dfun) ){
    Y_cand = Y0
    return(.plausibility_contour(Y_cand,s_aug,Y_aug,w_aug,d_aug,pred_fun,thetaHat,T_dfun)$p_y)
  } else {
    p_df   = .plausibility_contour(Y_cand=NULL,s_aug,Y_aug,w_aug,d_aug,pred_fun,thetaHat,T_dfun)
    Y_cand = p_df$Y_cand
    p_y    = p_df$p_y

    Y0    = sapply(Y0, function(x) min(max(Y_cand), x))
    p_Y0  = sapply(Y0, function(x) p_y[which(Y_cand >= x)[1]])
    return( p_Y0 )
  }


}
