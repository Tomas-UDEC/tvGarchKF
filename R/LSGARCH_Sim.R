#' @title Series simulation of a LSGARCH model.
#' @description Function to simulate a series of a long n LSGARCH model.
#' @param n Long to simulate the series
#' @param gamma Vector containing coefficents of c.
#' @param alpha Vector containing coefficents of alpha.
#' @param beta Vector containing coefficents of beta.
#' @param type Vector of function type for c, alpha and beta.
#' @param exponentes Vector for exponenets in NoLineal.
#' @param trig Type of trigonometric function.
#' @param arg Value of argument for the trigonometric function.
#' @return Vector with two componentes, where the first component is dependent variable and the second component is independent variable.
#' @importFrom stats na.exclude
#' @importFrom stats rnorm
#' @export
LSGARCH_Sim <- function(n, gamma, alpha, beta, type = c("polynomial", "NoLineal","trigonometric"),
                              exponentes = NULL, trig = NULL, arg = NULL){
  e<-rnorm(n)
  x<-rep(0, n)
  sigma<-rep(0,n)
  u<-(1:n)/n

  aux <- list(gamma=gamma, alpha=alpha, beta=beta)
  aux3 <- list()

  for(i in 1:length(type)){
    if(type[i] == "polynomial"){
      aux2 <- polynomial(n, coeficientes = unlist(aux[i]))
    }
    if(type[i] == "NoLineal"){
      aux2 <- NoLineal(n, coeficientes = unlist(aux[i]), exponentes = exponentes)
    }
    if(type[i] == "trigonometric"){
      aux2 <- trigonometric(n, coeficientes = unlist(aux[i]), trig = trig, arg = arg)
    }
    aux2 <- list(aux2)
    aux3 <- c(aux3,aux2)
  }

  gamma <- unlist(aux3[1])
  a <- unlist(aux3[2])
  b <- unlist(aux3[3])

  for(i in 2:n){
    sigma[i]<- gamma[i]+a[i]*x[i-1]**2+b[i]*sigma[i-1]
    x[i]<-e[i]*sqrt(sigma[i])
  }

  Y<- x^2
  for(i in 2:n){
    Y[i] <- gamma[i]*+(b[i]+a[i])*Y[i-1]-b[i]*e[i-1]+e[i]
  }
  Z<-cbind(x,Y)
  return(Z)
}
