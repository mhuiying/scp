#' Select the optimal tunning parameter for \code{\link{scp}}
#'
#' @param s an \eqn{n \times 2}{n x 2} \code{matrix} or a \code{data.frame} with two coordinates of \eqn{n} locations.
#' @param Y a vector with \eqn{n} values corresponding to \code{Y(s)}.
#' @param eta_cand a vector of candidate \code{eta} values. Defaults to \code{NULL}.
#' @param plot logical; if \code{TRUE}, \code{select_eta} plots the empirical interval score for the candidate \code{eta} values.
#'             if \code{FALSE}, no plot will be generated. Defaults to \code{TRUE}.
#'
#' @return The optimum \code{eta} based on minimizing empirical interval score.
#' @export
#'
#' @examples
#' #?sample_data
#' s  = sample_data$s
#' Y  = sample_data$Y
#'
#' # CAUTION: the following command may take one hour to run.
#' opt_eta = select_eta(s = s, Y = Y)
#'
#' use optimal eta to calculate prediction interval
#' s0  = c(0.5,0.5)
#' scp(s0 = s0, s = s, Y = Y, eta = opt_eta)

select_eta = function(s,Y,eta_cand = NULL, plot = TRUE){

  bins     = seq(0.01,0.2,0.01)
  thetaHat = get_theta(s,Y,dists=bins)

  valid_id = sample(1:nrow(s),100)
  intScore = c()

  min_d = min(dist(s))
  eta_cand = seq(3*min_d, 20*min_d, min_d)
  for(eta in eta_cand){
    PIs = scp(s[valid_id,],s=s,Y=Y,thetaHat=thetaHat,eta=eta)
    intScore = c(intScore, mean(interval_score(PIs[,1], PIs[,2], Y[valid_id])))
  }
  opt_eta = eta_cand[which.min(intScore)]

  if(plot){
    plot(eta_cand, intScore, type = "b", xlab = "eta", ylab = "interval score")
    abline(v = opt_eta, col = "red", lty=2)
  }

  return(opt_eta)

}
