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

# transform pred_fun to be suitable for plausibility fun
# trans_pred_fun = function(pred_fun,s_aug,Y_aug,Q=NULL){
#
#   if(!is.null(Q)){
#     Yhat = as.numeric(-Q%*%Y_aug/diag(Q))
#     sd   = as.numeric(sqrt(diag(Q)))
#   }else{
#     n    = length(Y_aug)
#     Yhat = rep(NA,n)
#     sd   = rep(1 ,n)
#     for(i in 1:n)
#       Yhat[i] = pred_fun(as.numeric(s_aug[i,]),s_aug[-i,],Y_aug[-i])
#   }
#
#   return(as.data.frame(Yhat=Yhat,sd=sd))
#
# }
