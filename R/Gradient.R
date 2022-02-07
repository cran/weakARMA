#' Computation the gradient of the residuals of an ARMA model
#'
#' @description Computes the gradient of the residuals of an ARMA model.
#'
#' @param ar Vector of \code{ar} coefficients.
#' @param ma Vector of \code{ma} coefficients.
#' @param y Univariate time series.
#' 
#'
#' 
#' @return A list containing:
#' \describe{
#'    \item{\code{der.eps}}{Matrix of the gradient.}
#'    \item{\code{esp}}{Vector of residuals.}
#' }
#' @export
#'
#' @examples 
#' est<-estimation(p = 1, q = 1, y = CAC40return.sq)
#' gradient(ar = est$ar, ma = est$ma, y = CAC40return.sq)
#' 
gradient <- function(ar = NULL, ma = NULL, y)
{
  
  grand = 1 / sqrt(.Machine$double.eps)
  n <- length(y)
  
  eps <- residual(ar = ar, ma = ma, y = y)
  C <- matrix(0, nrow = 0, ncol = n)
  if (is.null(ar))
  {
    p <- 0
    q <- length(ma)
    der.eps <- matrix(0, nrow = q, ncol = n)
    #Calcul de epsilon(t) = y(t)-a(1)y(t-1)-a(2)y(t-2)-...-a(p)y(t-p)-b(1)eps(t-1)-b(2)eps(t-2)-...-b(q)eps(t-q)
    
    
    #calcul de la derivee de eps / theta = -(y(t-1) ; y(t-2) ; ... ; y(t-p) ; (eps(t-1) ; eps(t-2) ; ... ; eps(t-q))-b(1)eps(t-1)-b(2)eps(t-2)-...-b(q)eps(t-q)
    for (t in 2:n){
      for (i in 1:min(q,(t-1))){
        der.eps[i, t] <- eps[t - i]
      }
    }
    
    C<-der.eps
    
    for(t in (q + 1):n) {
      for(i in 1:q){
        der.eps[ ,t] <- C[ ,t]+ma[i]*der.eps[ ,t - i]
      }
    }
  }
  else {
    if (is.null(ma)) {
      q <- 0
      p <- length(ar)
      der.eps <- matrix(0, nrow = p, ncol = n)
      
      #Calcul de epsilon(t) = y(t)-a(1)y(t-1)-a(2)y(t-2)-...-a(p)y(t-p)
      
      
      #calcul de la derivee de eps / theta = -(y(t-1) ; y(t-2) ; ... ; y(t-p))
      for (t in 2:n){
        for (i in 1:min(p,(t-1))){
          der.eps[i, t] <- -y[t - i]
        }
      }
    }
    else {
      p <- length(ar)
      q <- length(ma)
      der.eps.y <- matrix(0, nrow = p, ncol = n)
      der.eps.e <- matrix(0, nrow = q, ncol = n)
      der.eps <- matrix(0, nrow = (p + q), ncol = n)
      
      #Calcul de epsilon(t) = y(t)-a(1)y(t-1)-a(2)y(t-2)-...-a(p)y(t-p)-b(1)eps(t-1)-b(2)eps(t-2)-...-b(q)eps(t-q)
      
      
      #calcul de la derivee de eps / theta = -(y(t-1) ; y(t-2) ; ... ; y(t-p) ; (eps(t-1) ; eps(t-2) ; ... ; eps(t-q))-b(1)eps(t-1)-b(2)eps(t-2)-...-b(q)eps(t-q)
      for (t in 2:n){
        for (i in 1:min(q,(t-1))){
          der.eps.e[i, t] <- eps[t - i]
        }
        for (i in 1:min(p,(t-1))){
          der.eps.y[i, t] <- -y[t - i]
        }
      }
      
      C <- rbind(der.eps.y, der.eps.e)
 
      for(t in (q+1):n)
      {
        for(i in 1:(p+q)){
          der.eps[i,t] <- C[i,t]+sum(ma*der.eps[i,(t-1):(t-q)])
        }
      }
    }
  }
  
  return(list(gradient = der.eps, eps = eps))
}
