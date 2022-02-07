#'Computation of autocovariance and autocorrelation for an ARMA residuals.
#'
#' @description  Computes empirical autocovariances and autocorrelations
#'  functions for an ARMA process for only one given lag.
#'
#' @param ar Vector of AR coefficients. If \code{NULL},  it is a MA process.
#' @param ma Vector of MA coefficients. If \code{NULL},  it is a AR process.
#' @param y Univariate time series.
#' @param h Given lag to compute autocovariance and autocorrelation, with h an integer.
#' @param e Vector of residuals of the time series.  If \code{NULL}, the function will compute it.
#'
#' @return A list with :
#' \describe{
#'     \item{\code{autocov}}{Value of the autocovariance.}
#'     \item{\code{autocor}}{Value of the autocorrelation.}
#' }
#' 
#' @export
#'
#' @examples
#' param.estim <- estimation(p = 1,  q = 1, y = CAC40return.sq)
#' \donttest{acf.univ(ar = param.estim$ar, ma = param.estim$ma, y = CAC40return.sq,  h = 20)}
#' 
#'
#' @seealso \code{\link{acf.gamma_m}} for autocorrelation and autocovariance for all h lag.

acf.univ <- function (ar = NULL, ma = NULL, y, h, e = NULL)
{
  n <- length(y)
  if (is.null(e)){
    e <- residual(ar = ar, ma = ma,  y = y)
  }
  h <- floor(abs(h))
  centre <- e - mean(e)
  autocov <- sum(centre[1:(n - h)] * centre[(h + 1):n]) / n
  autocor <- autocov / (sum(centre[1:n] * centre[1:n]) / n)
  list(autocov = autocov, autocor = autocor)
}

#'Computation of autocovariance and autocorrelation for an ARMA residuals.
#'
#' @description  Computes empirical autocovariances and autocorrelations function
#'  for an ARMA process for lag max given.
#'
#' @param ar Vector of AR coefficients. If \code{NULL},  it is a MA process.
#' @param ma Vector of MA coefficients. If \code{NULL},  it is a AR process.
#' @param y Univariate time series.
#' @param h Computes autocovariances and autocorrelations from lag 1 to lag h with h an integer.
#' @param e Vector of residuals. If \code{NULL}, the function will compute it.
#'
#' @return A list with :
#' \describe{
#'     \item{\code{gamma_m}}{Vector of the autocovariances.}
#'     \item{\code{rho_m}}{Vector of the autocorrelations.}
#' }
#' 
#' 
#' @export
#'
#' @examples
#' param.estim <- estimation(p = 1,  q = 1, y = CAC40return.sq)
#' \donttest{acf.gamma_m(ar = param.estim$ar, ma = param.estim$ma, y = CAC40return.sq,  h = 20)}
#'
#' @seealso \code{\link{acf.univ}} for autocorrelation and autocovariance for only one given lag h.
#'
acf.gamma_m <- function (ar = NULL, ma = NULL, y, h, e = NULL)
{
  n <- length(y)
  if (is.null(e)){
    e <- residual(ar = ar, ma = ma,  y = y)
  }
  h <- abs(h)
  hmax <- min(n, h)
  autocov <- rep(0, hmax)
  autocor <- rep(0, hmax)
  centre <- e - mean(e)
  for (i in 1:hmax)
    autocov[i] <- sum(centre[1:(n - i)]*centre[(i + 1):n]) / n
  for (i in 1:hmax)
    autocor[i] <- autocov[i] / (sum(centre[1:n] * centre[1:n]) / n)
  list(gamma_m = autocov, rho_m = autocor)
}


#'Autocorrelogram
#'
#' @description Plots autocorrelogram for non linear process.
#'
#' @param ar Vector of AR coefficients. If \code{NULL}, we consider a MA process.
#' @param ma Vector of MA coefficients. If \code{NULL}, we consider an AR process.
#' @param y Univariate time series.
#' @param main Character string representing the title for the plot.
#' @param nlag Maximum lag at which to calculate the acf. If \code{NULL}, it is
#'  determinate by \eqn{nlag = min(10log(n))} where n is the number of
#'  observation.
#' @param conflevel Value of the confidence level, 5\% by default. 
#' @param z Zoom on the graph.
#' @param aff Specify the method between SN, M and both (see in Details).
#'
#' @importFrom stats qnorm acf
#' @importFrom graphics lines
#' 
#' @export
#' 
#' @examples 
#' est<-estimation(p = 1, q = 1, y = CAC40return.sq)
#' \donttest{nl.acf(ar = est$ar, ma = est$ma, y = CAC40return.sq, main = "Autocorrelation of an ARMA(1,1) 
#' residuals of the CAC40 return square", nlag = 20)}
#'
#' @note The only value available for the argument \code{conflevel} are 
#'  0.1, 0.05, 0.025, 0.01 or 0.005.
#'  
#' @details For the argument \code{aff} you have the choice between: 
#'  \code{SN}, \code{M} and \code{both}.
#'  \code{SN} prints the self-normalized method (see Boubacar Maïnassara and Saussereau) in green, 
#'  \code{M} prints the modified method introduced by Francq, Roy and Zakoïan (see also Boubacar Maïnassara) in red 
#'  and \code{both} prints both of the methods.
#'
#' @return An autocorrelogram with every autocorrelations from 1 to a lag max, and 
#' with methods you choose to print.
#'  
#' @references Boubacar Maïnassara, Y. 2011, Multivariate portmanteau test for structural {VARMA} models
#' with uncorrelated but non-independent error terms \emph{Journal of Statistical Planning and Inference},
#' vol. 141, no. 8, pp. 2961-2975.
#' @references Boubacar Maïnassara, Y.and Saussereau, B. 2018, Diagnostic checking in multivariate {ARMA} models with 
#' dependent errors using normalized residual autocorrelations ,
#'  \emph{Journal of the American Statistical Association}, vol. 113, no. 524, pp. 1813-1827.
#' @references Francq, C., Roy, R. and Zakoïan, J.M. 2005, Diagnostic Checking in ARMA
#'  Models with Uncorrelated Errors, \emph{Journal of the American Statistical
#'  Association}, vol. 100, no. 470, pp. 532-544.
#' @references Lobato, I.N. 2001, Testing that a dependant process is
#'  uncorrelated. J. Amer. Statist. Assos. 96, vol. 455, pp. 1066-1076.
#' 
#' 

nl.acf <- function(ar = NULL, ma = NULL, y, main = NULL, nlag = NULL, conflevel = 0.05, z=1.2, aff="both")
{
  if (conflevel!=0.1 & conflevel!=0.05 & conflevel!=0.025 & conflevel!=0.01 & conflevel!=0.005){
    stop("Choose a confidence level equal to 0.1, 0.05, 0.025, 0.01 or 0.005")
  }
  n <- length(y)
  if (is.null(nlag))
    nlag <- as.integer(min(10*log10(n), n - 1))

  coef.critique.SN <- 0
  coef.critique.ST <- 0
  auto <- acf.gamma_m(ar = ar, ma = ma, y = y, h = nlag) 
  rho_m <- auto$rho_m
  gamma_m <- auto$gamma_m

  if (is.null(ar) & is.null(ma)) res <- y
  #else res <- resultat.h(ar = ar, ma = ma, y = y, m = 1)$eps
  else res <- residual(ar = ar, ma = ma, y = y)
  
  var <- rep(0, nlag)
  var1 <- rep(0,nlag)
  for (h in 1:nlag) {
    mat <- resultat.h(ar = ar, ma = ma, y = y, m = h, eps=res)
    var[h] <- mat$Sigma_rho.h
    var1[h] <- mat$matC_m.h
  }

  band <- sqrt(var / n)
  band1 <- sqrt(var1 / n)
  
  if (aff!="M"){
    if (conflevel == 0.1)  {coef.critique.SN <- 28.31 }
    if (conflevel == 0.05) {coef.critique.SN <- 45.4 }
    if (conflevel == 0.025) {coef.critique.SN <- 66.13 }
    if (conflevel == 0.01) {coef.critique.SN <- 99.76 }
    if (conflevel == 0.005) {coef.critique.SN <- 128.1 }
  }
  if (aff!="SN"){
    coef.critique.ST<-abs(qnorm(conflevel/2))
  }

  
  minval <- z * min(rho_m, - sqrt(coef.critique.SN) * band1, - coef.critique.ST * band, - coef.critique.ST / sqrt(n))
  maxval <- z * max(rho_m, sqrt(coef.critique.SN) * band1, coef.critique.ST * band, coef.critique.ST / sqrt(n))

  acf(res, lag.max = nlag, xlab = 'Lag', ylab = 'ACF', ylim = c(minval, maxval), main = main, ci = (1-conflevel))
  
  if (aff == "M"){
    lines(c(1:nlag), -coef.critique.ST * band, lty = 1, col = 'red')
    lines(c(1:nlag), coef.critique.ST * band, lty = 1, col = 'red')
  }
  if (aff == "SN"){
    lines(c(1:nlag), -sqrt(coef.critique.SN) * band1, lty = 1, col = 'green')
    lines(c(1:nlag), sqrt(coef.critique.SN) * band1, lty = 1, col = 'green')
  }
  if (aff == "both"){
    lines(c(1:nlag), -coef.critique.ST * band, lty = 1, col = 'red')
    lines(c(1:nlag), coef.critique.ST * band, lty = 1, col = 'red')
    lines(c(1:nlag), -sqrt(coef.critique.SN) * band1, lty = 1, col = 'green')
    lines(c(1:nlag), sqrt(coef.critique.SN) * band1, lty = 1, col = 'green')
  }
  
 
}
