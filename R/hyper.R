#' Transform the hyperparameters
#'
#' @param s0
#' @param s
#' @param Y
#' @param global
#' @param eta
#' @param m
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
.hyper = function(s0,s,Y,global,eta,m){

  if(!global){
    dist = sqrt( (s0[1]-s[,1])^2 + (s0[2]-s[,2])^2 )
    these = (dist < 2*eta)
    dist_mat = apply(s[these,], 1, FUN = function(s0) sqrt( (s0[1]-s[,1])^2 + (s0[2]-s[,2])^2 ) )
    these = apply(dist_mat, 2, FUN = function(x) order(x)[1:15])
    these = unique(as.vector(these))
    if(is.null(m)) m = length(these)
  } else{
    m = length(Y)
    these = 1:m
  }

  s_aug = rbind(s0, s[these,])
  Y_aug = c(NA, Y[these])
  d_aug = as.matrix(dist(s_aug))
  if(eta==Inf){
    w = rep(1/(m+1),m+1)
  }else{
    w = d[1,-1]/eta
    w = c(1, exp(-0.5*w^2))
    w = w/sum(w)
  }

  sYw_aug = list(s_aug = s_aug, Y_aug = Y_aug, w_aug = w)
  return(sYw_aug)

}
