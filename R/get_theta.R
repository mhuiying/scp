#' Estimate theta via variogram fitting
#'
#' @param s spatial location
#' @param Y data
#' @param dists breakpoints for bins
#' @param plot_fitted will plot a emperical variogram if specified as TRUE
#'
#' @return Matern covariance parameters
#' @export
#'
#' @importFrom geoR variog
#' @importFrom dplyr left_join
#' @keywords internal
#'
get_theta = function(s,Y,dists=NULL){

  if( length(Y) != nrow(s) )
    stop( paste("the number of Y obs,", length(Y), ", does not match the number of locations,", nrow(s)) )

  if(length(Y) < 10^6){
    vario = variog(coords=s,data=Y,uvec=dists,messages=FALSE)
  } else {
    s1       = sort(unique(s[,1]))
    s2       = sort(unique(s[,2]))
    s_full   = expand.grid(x=s1,y=s2)
    sY       = as.data.frame(cbind(s, Y))
    sY_full  = left_join(s_full, sY, by = c("x" = "x", "y" = "y"))
    s        = sY_full[, 1:2]
    Y        = sY_full[, 3]

    vario = list( u=c(), v=c(), n=c() )
    s1    = sort(unique(s[,1]))
    n_s1  = length(s1)
    diff_s1 = s1[2]-s1[1]
    s2    = sort(unique(s[,2]))
    n_s2  = length(s2)
    Y_mat = matrix(Y, nrow = n_s1) # Y in matrix format

    for( NN in 1:10){
      D = c(Y_mat[1:(n_s1-NN),]-Y_mat[(NN+1):n_s1,], Y_mat[,1:(n_s2-NN)]-Y_mat[,(NN+1):n_s2])
      D = D[!is.na(D)]
      vario$n = c( vario$n, length(D) )
      vario$u = c( vario$u, diff_s1*NN )
      vario$v = c( vario$v, mean(D^2) )
    }
  }

  wsse = function(logtheta,u,v,w){
    theta = exp(logtheta)
    m     = theta[1]+theta[2]-mat_cov(u,exp(logtheta))
    out   = sum(w*(v-m)^2)
    return(out)}

  init  = c(0.5*var(Y,na.rm=TRUE),0.5*var(Y,na.rm=TRUE),median(vario$u),0.5)
  fit   = optim(log(init),wsse,u=vario$u,v=vario$v,w=vario$n)
  theta = exp(fit$par)
  names(theta) = c("Nugget","PartialSill","Range","Smoothness")

  return(theta)}
