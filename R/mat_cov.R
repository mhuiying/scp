#' Matern covariance function
#'
#' @param d distances
#' @param theta covariance parameters
#'
#' @return a covariance matrix
#' @export

mat_cov    <- function(d,theta){
  suppressMessages(library(geoR))
  theta[1]*(d==0) + theta[2]*matern(d,theta[3],theta[4])
}
