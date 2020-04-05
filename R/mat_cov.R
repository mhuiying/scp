#' Matern covariance function
#'
#' @param d a numeric distance, a vector of distances, or a distance matrix
#' @param theta Matern covariance parameters
#'
#' @return a numeric covariance, a vector of covariances, or a covariance matrix in the same size of \eqn{d}
#' @export
#'
#' @importFrom geoR matern

mat_cov = function(d,theta){
  theta[1]*(d==0) + theta[2]*matern(d,theta[3],theta[4])
}
