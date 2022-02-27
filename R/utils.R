#' Matern covariance function
#'
#' @param d a numeric distance, a vector of distances, or a distance matrix
#' @param theta Matern covariance parameters
#'
#' @return a numeric covariance, a vector of covariances, or a covariance matrix in the same size of \eqn{d}
#' @export
#'
#' @importFrom geoR matern
#'
#' @keywords internal
#'
mat_cov = function(d,theta){
  theta[1]*diag(nrow(d)) + theta[2]*matern(d,theta[3],theta[4])
}


#' Calculate interval score
#'
#' @param l interval lower bound
#' @param u interval upper bound
#' @param x data
#' @param alpha significance level. Defaults to 0.05.
#'
#' @return a value of interval score
#' @export
#'
#' @references \url{https://www.stat.washington.edu/raftery/Research/PDF/Gneiting2007jasa.pdf}, section 6.2
#' @keywords internal
#'
#' @examples
#' eps = 0.05
#' x = rnorm(1000)
#' u = qnorm(1-eps/2)
#' l = qnorm(eps/2)
#' interval_score(l, u, x, alpha = eps)
#'
interval_score = function(l, u, x, alpha = 0.05){
  width = u - l
  sl = 2 / alpha * (l - x) * (x < l)
  su = 2 / alpha * (x - u) * (x > u)
  width + sl + su
}
