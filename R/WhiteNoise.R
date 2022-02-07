#' Weak white noise
#' @description Simulates an uncorrelated but dependant noise process.
#' @param n Number of observations.
#' @param sigma Standard deviation.
#' @param k Integer \eqn{\neq 0} to prevent a zero denominator.
#' @param ninit Length of 'burn-in' period.
#' 
#' @importFrom stats rnorm
#' 
#' @export
#' 
#' @return Vector of size \code{n} containing  a nonlinear sequence \eqn{X_{i}} such as
#'   \eqn{X_i = \frac{Z_{i}}{|Z_{i+1}| + k}} , where \eqn{Z_{i}} is a sequence of iid
#'     random variables mean-zero random variable with variance \eqn{\sigma^2}.
#' @seealso \code{\link{wnPT}}, \code{\link{wnPT_SQ}}, \code{\link{simGARCH}}
#'
#' @references Romano, J. and Thombs, L. 1996, Inference for autocorrelation under weak assumptions,
#'  \emph{Journal of the American Statistical Association}, vol. 91, no. 434, pp. 590-600
#'
#' @examples
#' wnRT(100)
#' wnRT(100, sigma = 1)
wnRT <- function (n, sigma = 1, k = 1, ninit = 100)
{
  eps <- rep(0, (n + ninit))
  eta <- rnorm((n + ninit + 1), sd = sigma)
  for (t in 1:(n + ninit)) {
    eps[t] <- eta[t + 1] / (abs(eta[t]) + k)}
  return(eps[(ninit + 1): (n + ninit)])
}


#' Weak white noise
#' @description Simulates an uncorrelated but dependant noise process.
#' @param n Number of observations.
#' @param sigma Standard deviation.
#' @param k Integer corresponding to the number of past observation will be used.
#' @param ninit Length of 'burn-in' period.
#' 
#' @importFrom stats rnorm
#' 
#' @export
#' 
#' @return Vector of size \code{n} containing  a nonlinear sequence \eqn{X_{i}} such as
#'   \eqn{X_{i} = Z_{i}Z_{i-1}...Z_{i-k}} , where \eqn{Z_{i}} is a sequence of iid
#'     random variables mean-zero random variable with variance \eqn{\sigma^2}.
#' @seealso \code{\link{wnRT}}, \code{\link{wnPT_SQ}}, \code{\link{simGARCH}}
#'
#' @references Romano, J. and Thombs, L. 1996, Inference for autocorrelation under weak assumptions,
#'  \emph{Journal of the American Statistical Association}, vol. 91, no. 434, pp. 590-600
#'
#' @examples
#' wnPT(100)
#' wnPT(100, sigma = 1, k = 1)
#' wnPT(100, k = 0) #strong noise
wnPT <- function (n, sigma = 1, k = 1, ninit = 100)
{
  eps <- rep(0, (n + ninit))
  eta <- rnorm((n + ninit) , sd = sigma)
  for(t in (k + 1):(n + ninit))
    eps[t] <- prod(eta[t:(t - k)])
  return(eps[(ninit + 1) : (n + ninit)])
}

#' Weak white noise
#' @description Simulates an uncorrelated but dependant noise process.
#' @param n Number of observations.
#' @param sigma Standard deviation.
#' @param k Integer corresponding to the number of past observation will be used.
#' @param ninit Length of 'burn-in' period.
#' 
#' @importFrom stats rnorm
#' @export
#' @return Vector of size \code{n} containing  a nonlinear sequence \eqn{X_{i}} such as
#'   \eqn{X_{i} = Z^{2}_iZ_{i-1}...Z_{i-k}} , where \eqn{Z_{i}} is a sequence of iid
#'     random variables mean-zero random variable with variance \eqn{\sigma^2}.
#' @seealso \code{\link{wnRT}}, \code{\link{wnPT}}, \code{\link{simGARCH}}
#'
#' @references Romano, J. and Thombs, L. 1996, Inference for autocorrelation under weak assumptions,
#'  \emph{Journal of the American Statistical Association}, vol. 91, no. 434, pp. 590-600
#'
#' @examples
#' wnPT_SQ(100)
#' wnPT_SQ(100, sigma = 1, k = 1)
wnPT_SQ <- function (n, sigma = 1, k = 1, ninit=100)
{
  eps <- rep(0, (n + ninit))
  eta <- rnorm(n + ninit, sd = sigma)
  for(t in (k + 1):(n + ninit))
    eps[t] <- eta[t]*prod(eta[t:(t - k)])
  return(eps[(ninit + 1) : (n + ninit)])
}


#' GARCH process
#' @description Simulates a GARCH process which is an example of a weak white noise.
#' 
#' @param n Number of observations.
#' @param c Positive number.
#' @param A Vector of ARCH coefficients >=0.
#' @param B Vector of GARCH coefficients >=0. If \code{NULL}, the
#'   simulation is a ARCH process.
#' @param ninit Length of 'burn-in' period.
#' 
#' 
#' @importFrom stats rnorm
#' @export
#' 
#' @return Vector of size \code{n} containing  a nonlinear sequence \eqn{\epsilon_t} such as
#'   \deqn{\epsilon_{t} = H_{t}^{1 / 2}  \eta_{t}} where \deqn{H_{t} =  c +
#'   a_{1}\epsilon_{t - 1}^ {2}+...+a_{q}\epsilon_{t - q} ^{2} + b_{1}H_{t-1}+...+ b_{p}H_{t-p}}
#'
#' @references Francq C. and Zakoïan J.M., 2010, \emph{GARCH models: structure, statistical inference and financial applications}
#'  
#'  
#' @seealso \code{\link{wnRT}}, \code{\link{wnPT}}, \code{\link{wnPT_SQ}}
#'
#' @examples
#' simGARCH(100, c = 1, A = 0.25)
#' simGARCH(100, c = 1, A = 0.1,  B = 0.88)


simGARCH <- function(n, c, A, B = NULL, ninit = 100)
{
  q <- length(A)
  p <- length(B)

  if (missing(B)) {B <- 0 ; p <- 0}

  eps <- rep(0 , n + ninit)
  eta <- rnorm(n + ninit, sd = 1)
  H <- rep(0, n + ninit)
  H[1] <- c
  eps[1] <- H[1] * (eta[1])
  
  
  for (t in 2:(n + ninit)){
    H[t] <- c + sum(A[1:min((t-1),q)] * (eps[(t-1):max(1,(t-q))])^2) + sum(B[1:min((t-1),p)] * H[(t-1):max(1,(t-p))])
    eps[t] <- sqrt(H[t]) * eta[t]
  }

  eps <- eps[(1 + ninit):(n + ninit)]
  return(eps)
}
