#' Transform the input parameters
#'
#' @description An internal function to transform the input parameters
#'
#' @param s0 prediction location, a numeric vector with \code{length = 2},
#'                             or a \code{matrix} with 1 row and 2 cols,
#'                             or a data.frame with 1 row and 2 cordinates.
#' @param s an \eqn{n \times 2}{n x 2} \code{matrix} or \code{data.frame} with two coordinates of \eqn{n} locations.
#' @param Y a vector with \eqn{n} values corresponding to \code{Y(s)}.
#' @param global logical; if \code{TRUE} , \code{scp} function returns the result of global spatial conformal prediction (GSCP);
#' if \code{FALSE}, \code{scp} function returns the result of local spatial conformal prediction (LSCP)
#' and users need to specify \code{eta}. Defaults to \code{TRUE}.
#' @param eta kernel bandwidth for weight schema, a positve scalar with smaller value meaning more localized procedure.
#' Defauls to \code{Inf}, which puts equal weight on surrounding \eqn{m} points.
#' @param m an postive integer representing the number of nearest locations to use for prediction.
#' Default depands on \code{eta}.
#' @param dfun non-conformity measure with four options.
#'             In which, \code{"residual2"} (default) represents squared residual,
#'             \code{"std_residual2"} represents standardized squared residual,
#'             \code{"abs_residual"} represents absolute residual,
#'             and \code{"std_abs_residual"} represents standardized absolute residual.
#'
#' @export
#' @keywords internal

.prime = function(s0,s,Y,global,eta,m,dfun){

  # idx = which( s[,1] == as.numeric(s0[1]) & s[,2] == as.numeric(s0[2]) )
  # if(length(idx) > 0){s = s[-idx,]; Y = Y[-idx]}

  these = .deter_these(s0,s,Y,global,eta,m)
  s_aug = rbind(s0, s[these,])
  Y_aug = c(NA, Y[these])
  d_aug = as.matrix(dist(s_aug))

  M = length(these)
  if(eta==Inf){
    w = rep(1/(M+1),M+1)
  }else{
    w = d_aug[1,-1]/eta
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

#' Determine the surrounding "these" points for prediction
#'
#' @param s0 prediction location, a numeric vector with \code{length = 2}.
#' @param s an \eqn{n \times 2}{n x 2} \code{matrix} or \code{data.frame} with two coordinates of \eqn{n} locations.
#' @param Y a vector with \eqn{n} values corresponding to \code{Y(s)}.
#' @param global logical; if \code{TRUE} , \code{scp} function returns the result of global spatial conformal prediction (GSCP);
#' if \code{FALSE}, \code{scp} function returns the result of local spatial conformal prediction (LSCP)
#' and users need to specify \code{eta}. Defaults to \code{TRUE}.
#' @param eta kernel bandwidth for weight schema, a positve scalar with smaller value meaning more localized procedure.
#' Defauls to \code{Inf}, which puts equal weight on surrounding \eqn{m} points.
#' @param m an postive integer representing the number of nearest locations to use for prediction.
#' Default depands on \code{eta}.
#'
#' @return a vector of location indexes for prediction
#' @export
#'
#' @keywords internal
.deter_these = function(s0,s,Y,global,eta,m){

  if( eta <= 0 ){
    stop("Please provide a positive eta")
  }  else if( eta < Inf ){
    global = FALSE
  }

  if( !is.null(m) ){
    if( (m <= 0 | m%%1!=0) ){
      stop("Please provide a positive integer for m")
    }
    else if( m < length(Y) ){
      global = FALSE
    }
  }


  if(global){ # GSCP

    these = 1:length(Y)

  } else {    # LSCP

    dist = sqrt( (as.numeric(s0[1])-s[,1])^2 + (as.numeric(s0[2])-s[,2])^2 )
    if(is.null(m)){
      these = (dist < 2*eta)
      if( sum(these) > 5 ){
        dist_mat = apply(s[these,], 1, FUN = function(s0) sqrt( (as.numeric(s0[1])-s[,1])^2 + (as.numeric(s0[2])-s[,2])^2 ) )
        these = apply(dist_mat, 2, FUN = function(x) order(x)[1:15])
        these = unique(as.vector(these))
      } else {
        warning(paste("eta is too small. Try eta >", sort(dist)[5]/2))
        these = order(dist)[1:30]
      }
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
#' @param dfun non-conformity measure with four options.
#'             In which, \code{"residual2"} (default) represents squared residual,
#'             \code{"std_residual2"} represents standardized squared residual,
#'             \code{"abs_residual"} represents absolute residual,
#'             and \code{"std_abs_residual"} represents standardized absolute residual.
#'
#' @return
#' @export
#' @keywords internal

trans_dfun = function(dfun){

  if(dfun == "std_residual2"){
    dfun1 = function(y,yhat,sd) (y-yhat)^2/sd^2
    std = TRUE
  }else if(dfun == "residual2"){
    dfun1 = function(y,yhat,sd=1) (y-yhat)^2/sd^2
    std = FALSE
  }
  return(list(fun = dfun1, std = std, residual2 = grepl("residual2",dfun)))

}

#' Generate candidate Y values
#'
#' @param pred_fun spatial prediction function with inputs being \eqn{s0, s, Y} and ouputs being predicted \code{Y(s0)}
#' (and its standard error). Defaults to \code{\link{krige_pred}} representing Kriging prediction.
#' @param dfun non-conformity measure with four options.
#'             In which, \code{"residual2"} (default) represents squared residual,
#'             \code{"std_residual2"} represents standardized squared residual,
#'             \code{"abs_residual"} represents absolute residual,
#'             and \code{"std_abs_residual"} represents standardized absolute residual.
#' @param precision a positive scalar represents how dense the candidates for \code{Y(s)} are. Defaults to \code{NULL}.
#'
#' @return a vector of candidate Y values
#' @export
#'
#' @keywords internal
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
