#' @title Sample data to demonstrate function usage in \code{scp}
#'
#' @description A list containing locations, \code{s},
#' and corresponding observations, \code{Y}, to demostrating \code{\link{scp}},
#' \code{\link{plausibility}}, and \code{\link{plausibility_contour}} usage
#'
#' @format A \code{list} with 2 elements, which are:
#' \describe{
#' \item{s}{a \code{data.frame} with two coordinates of 1681 grid locations in \eqn{[0,1]^2}.}
#' \item{Y}{a vector with \eqn{Y(s)} observations corresponding to \code{s}, where \eqn{Y(s) = X(s)^3 + E(s)},
#' \eqn{X(s)} is a stationary Gaussian process process with a Matern covariance, and \eqn{E(s)} is a white noise process.}
#' }
"sample_data"
