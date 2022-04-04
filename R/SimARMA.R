#'Simulation of ARMA(p,q) model.
#' @description Simulates an ARMA, AR or MA process according to the arguments
#'  given.
#'
#' @param n Number of observations.
#' @param ar Vector of AR coefficients. If \code{NULL},  the simulation is a MA process.
#' @param ma Vector of MA coefficients. If \code{NULL},  the simulation is a AR process.
#' @param sigma Standard deviation.
#' @param eta Vector of white noise sequence. Allows the user to use his own
#'  white noise.
#' @param method Defines the kind of noise used for the simulation. By default,
#'  the noise used is strong. See 'Details'.
#' @param k Integer used in the creation of the noise. See 'Details'.
#' @param mu Integer for the mean of the series.
#' @param ... Arguments needed to simulate GARCH noise. See 'Details'.
#'
#' @details ARMA model is of the following form : \deqn{ X_{t}-\mu = e_{t} + a_{1} (X_{t-1}-\mu)
#' + a_{2} (X_{t-2}-\mu) + ... + a_{p} (X_{t-p}-\mu) - b_1 e_{t-1} - b_2 e_{t-2} - ... - b_{q}  e_{t-q}}
#'  where \eqn{e_t} is a sequence  of uncorrelated random variables with zero
#'  mean and common variance \eqn{\sigma^{2} > 0} . \eqn{ar = (a_{1}, a_{2}, ..., a_{p})} are
#'  autoregressive coefficients and \eqn{ma = (b_{1}, b_{2}, ... , b_{q})} are moving
#'  average coefficients. Characteristic polynomials of ar and ma must
#'  constitute a stationary process.
#'
#'  Method "\code{strong}" realise a simulation with gaussian white noise.
#'
#'  Method "\code{product}", "\code{ratio}" and "\code{product.square}"
#'  realise a simulation with a weak white noise. These methods employ
#'  respectively the functions  \code{\link{wnPT}}, \code{\link{wnRT}} and
#'  \code{\link{wnPT_SQ}} to simulate nonlinear ARMA model. So, the
#'  paramater \code{k} is an argument of these functions. See \code{\link{wnPT}}, \code{\link{wnRT}}
#'  or \code{\link{wnPT_SQ}}.
#'
#'  Method "\code{GARCH}" gives an ARMA process with a GARCH noise. See
#'  \code{\link{simGARCH}}.
#'
#' @return Returns a vector containing the \code{n} simulated observations of the
#'  time series.
#'
#' @importFrom stats rnorm
#' 
#' @export
#'
#' @examples
#' y <- sim.ARMA(n = 100, ar = 0.95, ma = -0.6, method = "strong" )
#' y2 <- sim.ARMA(n = 100, ar = 0.95, ma = -0.6, method = "ratio")
#' y3 <- sim.ARMA(n = 100,  ar = 0.95, ma = -0.6, method = "GARCH", c = 1, A = 0.1, B = 0.88)
#' y4 <- sim.ARMA(n = 100, ar = 0.95, ma = -0.6, method = "product")
#' y5 <- sim.ARMA(n = 100, ar = 0.95, ma = -0.6, method = "product.square")
#'
#' @references Francq, C. and Zakoïan, J.M. 1998, Estimating linear representations
#'  of nonlinear processes, \emph{Journal of Statistical Planning and
#'  Inference}, vol. 68, no. 1, pp. 145-165
#'
#' @seealso \code{\link[stats]{arima.sim}}
#'
sim.ARMA <- function (n, ar = NULL, ma = NULL, sigma = 1, eta = NULL, method = "strong", k = 1, mu=0, ...)
{
  p <- length(ar)
  q <- length(ma)
  
  if (is.null(eta)) {
    switch(method,
           "strong" =
             {eta <- rnorm(n, sd = sqrt(sigma))},
           "product" =
             {eta <- wnPT(n , sigma, k = k)},
           "product.square" =
             {eta <- wnPT_SQ(n, sigma, k = k)},
           "ratio" =
             {eta <- wnRT(n , sigma, k = k)},
           "GARCH" =
             {eta <- simGARCH(n = n , ...)}
    )
  }
  x <- eta
  xbis<-x-mu
  
  if (!is.null(ma)) if (min(abs(polyroot(c(1, -1*ma)))) < 1) warning ("Polynomial MA not invertible")
  if (!is.null(ar)) if (min(abs(polyroot(c(1, -1*ar)))) < 1) warning ("Polynomial AR not invertible")
  
  if (is.null(ma)) {ma <- 0 ; q <- 0}
  if (is.null(ar)) {ar <- 0 ; p <- 0}
  
  for (t in 2:n){
    xbis[t] <- eta[t] + sum(ar[1:min((t-1),p)]*xbis[(t - 1):max(1,(t - p))]) - sum(ma[1:min((t-1),q)]*eta[(t - 1):max(1,(t - q))]) #coherence avec article -> -sum(ar*x[(t-1):(t-p)]
  } 
  x<-xbis+mu
  return(x)
  
}
