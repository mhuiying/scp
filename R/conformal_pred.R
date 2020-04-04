#' Spatial Conformal Prediction Intervals
#'
#' @description This function provides the conformal prediction interval for spatial location \code{s0} given obserations \code{s, Y}.
#'
#' @param s0 prediction location
#' @param s an \eqn{n \times 2}{n x 2} matrix or data-frame with the 2-D coordinates of the \eqn{n} data locations.
#' @param Y a vector with \eqn{n} data values.
#' @param theta spatial covariance parameters as in \code{\link{mat_cov}}
#' @param r numerical value of the kernel bandwidth for the weight schema in conformal prediction. The defaul is \eqn{Inf} meaning equal weight on surrounding \eqn{m} points.
#' @param eps significance level. The default is 0.05.
#'
#'
#' @return A vector of lower and upper bounds of the conformal prediction interval.
#' @export
#'
#' @examples
#' N = 40; n = N^2
#' S = seq(0,1,length=40)
#' s <- expand.grid(S,S)
#' d <- as.matrix(dist(s))
#'
#' theta        = c(0,3,0.1,0.7)
#' names(theta) = c("Nugget","PartialSill","Range","Smoothness")
#' C <- mat_cov(d,theta)
#' X <- t(chol(C))%*%rnorm(n)
#' Y <- X^3 + rnorm(n)
#'
#' # Estimate spatial covariance parameters
#' bins     <- seq(0.01,0.2,0.01)
#' thetaHat <- get_theta(s,Y,dists=bins,plot_fitted=FALSE)
#' Q        <- solve(mat_cov(d,thetaHat))
#'
#' s0       <- c(0.5, 0.5)
#' conformal_pred(s0,s,Y,thetaHat,m=100,eta=0.1,eps=eps)

conformal_pred = function(s0,s,Y,theta,eta=Inf,eps=0.05){

  dist = sqrt( (s0[1]-s[,1])^2 + (s0[2]-s[,2])^2 )
  s1 = unique(s[,1])
  these <- (dist < 2*eta + 3/length(s1))
  m = sum(these)
  s     <- rbind(s0,s[these,])
  Y     <- c(NA,Y[these])
  d     <- as.matrix(dist(s))
  if(eta==Inf){
    w <- rep(1/m,m)
  }else{
    w <- d[1,-1]/eta
    w <- exp(-0.5*w^2)
    w <- w/sum(w)
  }

  Q      <- solve(mat_cov(d,theta))
  UVW    <- compute_UVW(Q,Y)
  U      <- UVW$U[-1]
  V      <- UVW$V[-1]
  W      <- UVW$W[-1]

  low    <- (-V + sqrt(V^2 - 4*U*W))/(2*W)
  upp    <- (-V - sqrt(V^2 - 4*U*W))/(2*W)
  lu_w   <- c(w,w)
  lu_lab <- c(rep(1, length(low)), rep(-1, length(upp)))
  lu_w   <- lu_w[order(c(low, upp))]
  lu_lab <- lu_lab[order(c(low, upp))]

  y0     <- sort(c(low, upp))
  p_y    <- cumsum(lu_w*lu_lab)

  gamma  <- range(y0[which(p_y >= floor((m+1) * eps) / (m+1))])
  return(gamma)}
