#' Function optim will minimize
#'
#' @param x One point in \eqn{\rm I\!R^{(p+q)}}.
#' @param dim.ar Length of AR vector.
#' @param dim.ma Length of MA vector.
#' @param y Vector of a time series.
#'
#' @return 
#' \describe{
#'    \item{\code{ms}}{Mean square at the point \code{x}.}
#' }
#' 
#' 
#' @description Computes the mean square of the time series at the point \code{x}, will be minimize with the \code{\link[stats]{optim}} function in our function \code{\link{estimation}}.
#' 
#' @export
#'

meansq <- function(x, dim.ar = NULL, dim.ma = NULL, y) 
{
  n <- length(y)
  e <- y
  if (is.null(dim.ar) || dim.ar==0) {
    ma <-x 
    ar <- 0 
    dim.ar <- 0
    for (t in  2:n){
      e[t] <- y[t]+ sum(ma[1:min((t-1),dim.ma)]*e[(t-1):max(1,(t - dim.ma))])
    }
  } else {
      if (is.null(dim.ma) || dim.ma==0) {
        ar <- x 
        ma <- 0 
        dim.ma <- 0
        for (t in  2:n){
          e[t] <- y[t]- sum(ar[1:min((t-1),dim.ar)]*y[(t-1):max(1,(t - dim.ar))]) 
        }
      } else {
          ar <- x[1:dim.ar]
          ma <- x[(dim.ar + 1):length(x)]
          for (t in  2:n){
            e[t] <- y[t]- sum(ar[1:min((t-1),dim.ar)]*y[(t-1):max(1,(t - dim.ar))]) + sum(ma[1:min((t-1),dim.ma)]*e[(t-1):max(1,(t - dim.ma))])
          }
        }
    }
  ms <- mean(e^2)
  return(ms)
}


#' Parameters estimation of a time series.
#' 
#' @description Estimates the parameters of a time series for given orders \code{p} and \code{q}
#'
#' @param p Order of AR, if \code{NULL}, MA is computed.
#' @param q Order of MA, if \code{NULL}, AR is computed.
#' @param y Univariate time series.
#' @param meanparam Logical argument if the mean parameter has to be computed or not. If FALSE \eqn{\mu} is not computed.
#' 
#' @importFrom stats optim
#'
#' @return List of estimate coefficients:
#' \describe{
#'     \item{\code{mu}}{Mean parameter}.
#'     \item{\code{ar}}{Vector of AR coefficients with length is equal to \code{p}.}
#'     \item{\code{ma}}{Vector of MA coefficients with length is equal to \code{q}.}
#'     \item{\code{sigma.carre}}{Mean square residuals.}
#' }
#' @export
#' 
#' @details This function uses the algorithm BFGS in the function optim to minimize our objective function \code{\link{meansq}}.
#'
#' @references Francq, C. and Zakoïan, J. 1998, Estimating linear representations of nonlinear processes
#' \emph{Journal of Statistical Planning and Inference}, vol. 68, no. 1, pp. 145-165.
#'
#' @examples 
#' y<-sim.ARMA(1000,ar = c(0.9,-0.3), ma = 0.2, method = "product")
#' estimation(p = 2, q = 1, y = y)
#' 
#' estimation(p = 1, q = 1, y = CAC40return.sq, meanparam = TRUE)
#' 
estimation <- function(p = NULL, q = NULL, y, meanparam = FALSE)
{
  mu<-0
  if (meanparam == TRUE){
    mu<-mean(y)
    y<-y-mu
  }
  
  if (is.null(p)) {p1<-0;q1<-q}
  else {
    if (is.null(q)) {q1<-0;p1<-p}
    else {p1<-p;q1<-q}
  } #si p et q ne peuvent pas etre egaux à 0 en meme temps
  
  para<-numeric(p1+q1)
  res <- optim(par = para, fn = meansq, dim.ar = p1, dim.ma = q1, y = y, method = "BFGS",control = list(trace=0))
  
  if(is.null(p) || p==0) {
    return(list(mu=mu,ma=res$par,sigma.carre=res$value))
  }
  else {
    if (is.null(q) || q==0){
      return(list(mu=mu,ar=res$par,sigma.carre=res$value))
    }
    else {
      return(list(mu=mu,ar=res$par[1:p],ma=res$par[(p+1):length(res$par)],sigma.carre=res$value))
    }
  }
} 


