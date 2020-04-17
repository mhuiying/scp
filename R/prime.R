#' Transform the hyperparameters
#'
#' @param s0
#' @param s
#' @param Y
#' @param global
#' @param eta
#' @param m
#' @param dfun
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
.prime = function(s0,s,Y,global,eta,m,dfun){

  idx = which(s[,1]==s0[1] & s[,2]==s0[2])
  if(length(idx) > 0){s = s[-idx,]; Y = Y[-idx]}

  these = .deter_these(s0,s,Y,global,eta,m)
  s_aug = rbind(s0, s[these,])
  Y_aug = c(NA, Y[these])
  d_aug = as.matrix(dist(s_aug))

  M = length(these)
  if(eta==Inf){
    w = rep(1/(M+1),M+1)
  }else{
    w = d[1,-1]/eta
    w = c(1, exp(-0.5*w^2))
    w = w/sum(w)
  }

  s_aug  <<- s_aug
  Y_aug  <<- Y_aug
  w_aug  <<- w
  d_aug  <<- d_aug
  M      <<- M

  T_dfun <<- trans_dfun(dfun)

}

.deter_these = function(s0,s,Y,global,eta,m){

  if( eta <= 0 )
    stop("Please provide a positive eta")
  else if( eta < Inf )
    global = FALSE

  if( !is.null(m) )
    if( (m <= 0 | m%%1!=0) )
      stop("Please provide a positive integer for m")
    else if( m < length(Y) )
      global = FALSE

  if(global){ # GSCP

    these = 1:length(Y)

  } else {    # LSCP

    dist = sqrt( (s0[1]-s[,1])^2 + (s0[2]-s[,2])^2 )
    if(is.null(m)){
      these = (dist < 2*eta)
      dist_mat = apply(s[these,], 1, FUN = function(s0) sqrt( (s0[1]-s[,1])^2 + (s0[2]-s[,2])^2 ) )
      these = apply(dist_mat, 2, FUN = function(x) order(x)[1:15])
      these = unique(as.vector(these))
    } else {
      these = order(dist)[1:m]
    }

  }

  if( !is.null(m) )
    if( m < 30)
      warning('m is smaller than 30. A less localized scp is recommended. Please either increase eta or m.')
    else if( m > 5000)
      warning('m is larger than 5000. A more localized scp is recommended. Please either decrease eta or m.')

  return(these)
}

#' Transform non-conformity measure
#'
#' @param dfun
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
trans_dfun = function(dfun){

  if(dfun == "std_residual2"){
    dfun1 = function(y,yhat,sd) (y-yhat)^2/sd^2
    std = TRUE
  }else if(dfun == "residual2"){
    dfun1 = function(y,yhat,sd=1) (y-yhat)^2/sd^2
    std = FALSE
  }else if(dfun == "abs_residual"){
    dfun1 = function(y,yhat,sd=1) abs(y-yhat)/sd
    std = FALSE
  }else if(dfun == "std_abs_residual"){
    dfun1 = function(y,yhat,sd) abs(y-yhat)/sd
    std = TRUE
  }else if(!is.function(dfun)){
    stop("Please provide a valid non-conformity measure.")
  }

  return(list(fun = dfun1, std = std, residual2 = grepl("residual2",dfun)))
}

.generate_Y_cand = function(pred_fun, dfun, precision){

  if( !identical(pred_fun, krige_pred) | !grepl("residual2", dfun) ){
    if(is.null(precision))
      Y_cand = seq(min(Y),max(Y),length.out = 100)
    else
      Y_cand = seq(min(Y),max(Y),precision)
  } else {
    Y_cand = NULL
  }
  return(Y_cand)

}
