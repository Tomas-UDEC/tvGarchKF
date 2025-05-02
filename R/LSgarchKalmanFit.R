#' @title Fit LSGARCH parameters with Kalman Filter.
#' @description Fit procedure for coefficents of c, \eqn{\alpha} and \eqn{\beta}.
#' @param series Time series.
#' @param c Vector containing coefficents of c.
#' @param alpha Vector containing coefficents of alpha.
#' @param beta Vector containing coefficents of beta.
#' @param type Vector of function type for c, alpha and beta.
#' @param exponentes Vector for exponenets in NoLineal.
#' @param trig Type of trigonometric function.
#' @param arg Value of argument for the trigonometric function.
#' @param predict Value for time to generate predict.
#' @param trace.log Variable to print names of coefficients.
#' @return Return fit values of omega, alpha and beta
#' @importFrom stats optim
#' @useDynLib LSGARCH, .registration = TRUE
#' @export
LSgarchKalmanFit <- function(series, c, alpha, beta, type=c("polynomial","NoLineal","trigonometric"), exponentes, trig, arg, predict = 0, trace.log=F){
  checkInput(type)
  pars = c(c, alpha, beta)
  if(!missing(predict)){
    predict <- predict
  }else{
    predict <- 0
  }
  if(type[1]!="polynomial" | type[2]!="polynomial" | type[3]!="polynomial"){
    aux <- optim(par = pars, fn = LSgarchKalmanLoglike, series = series, c = c, alpha = alpha, beta = beta, type = type, exponentes = exponentes, trig = trig, arg = arg, method = c("BFGS"))
  }else{
    aux <- optim(par = pars, fn = LSgarchKalmanLoglike, series = series, c = c, alpha = alpha, beta = beta, type = type, exponentes = exponentes, trig = trig, arg = arg, method = c("CG"))
  }
  par <- round(aux$par,5)
  if(trace.log == T){
    for(i in 1:max(length(c),length(alpha),length(beta))){
      if(i<=length(c)){
        names(pars)[i] <- paste0("c",i)
      }
      if(i<=length(alpha)){
        names(pars)[i+length(c)] <- paste0("a",i)
      }
      if(i<=length(beta)){
        names(pars)[i+length(c)+length(alpha)] <- paste0("b",i)
      }
    }
  }
  suppressWarnings(return(par))
}

