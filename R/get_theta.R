#' Estimate theta via variogram fitting
#'
#' @param s spatial location
#' @param Y data
#' @param dists breakpoints for bins
#' @param max_range maximum range
#'
#' @return Matern covariance parameters
#' @export
#'
#' @importFrom dplyr left_join
#' @keywords internal
#'
get_theta <- function(s, Y, dists = NULL, max_range=NULL) {

  if(is.null(max_range)){max_range=max(dist(s))}
  vario = variog(coords = s, data = Y, uvec = dists, messages = FALSE)

  wsse  = function(logtheta, u, v, w, max_range) {
    theta = exp(logtheta)
    theta[4] = max_range*pnorm(logtheta[4])
    m     = theta[1] + theta[2] - theta[2]*matern(u,theta[3],theta[4])
    out   = sum(w * (v - m)^2)
    return(out)}

  init  = c(0.5 * var(Y, na.rm = TRUE), 0.5 * var(Y, na.rm = TRUE),
            quantile(vario$u,0.5), 0.5)

  fit   = optim(log(init), wsse, u = vario$u, v = vario$v,
                w = vario$n, max_range=max_range)

  theta = exp(fit$par)

  theta[4] = max_range*pnorm(fit$par[4])

  names(theta) = c("Nugget", "PartialSill", "Range", "Smoothness")

  return(theta)}
