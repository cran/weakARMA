#' Computation of Fisher information  matrice 
#' 
#' @description Computes matrices of Fisher information like \eqn{I}, \eqn{J}.
#' 
#'
#' @param ar Vector of AR coefficients. If \code{NULL},  the simulation is a MA process.
#' @param ma Vector of MA coefficients. If \code{NULL},  the simulation is a AR process.
#' @param y Univariate time series.
#'
#' @return A list of matrix containing:
#' \describe{
#'    \item{\code{I}}{Matrix \code{I} computed in function \code{\link{matXi}}.}
#'    \item{\code{J}}{Matrix \code{J} computed as \eqn{\frac{2}{n} H(e) H(e)^t } where \eqn{e} is the residuals vector.}
#'    \item{\code{J.inv}}{Inverse of the matrix \code{J}.}
#'    \item{\code{matOmega}}{Matrix variance-covariance in the weak case computed as \eqn{J^{-1}IJ^{-1}}.}
#'    \item{\code{matvar.strong}}{Matrix variance-covariance in the strong case computed as 
#'          \eqn{2\sigma^2J^{-1}}.}
#'    \item{\code{standard.dev.Omega}}{Standard deviation of the matrix \code{matOmega}.}
#'    \item{\code{standard.dev.strong}}{Standard deviation of the matrix \code{matvar.strong}.}
#'    \item{\code{sig2}}{Innovation variance estimate.}
#' }
#' @import vars
#' 
#' @export
#'
#' @examples 
#' y <- sim.ARMA(n = 1000, ar = c(0.95,-0.8), ma = -0.6)
#' \donttest{est<-estimation(p = 2, q = 1, y = y)}
#' \donttest{omega(ar = est$ar, ma = est$ma, y = y)}
#' 
#' estCAC<-estimation(p = 1, q = 1, y = CAC40return.sq, meanparam = TRUE)
#' \donttest{omega(ar = estCAC$ar, ma = estCAC$ma, y = CAC40return.sq)}

omega <- function(ar = NULL, ma = NULL, y)
{
  
  grand = 1 / sqrt(.Machine$double.eps)
  n <- length(y)
  
  if (is.null(ma)) {p <- length(ar) ; q <- 0}
  else {
    if (is.null(ar)) {q <- length(ma) ; p <- 0}
    else {p <- length(ar) ; q <- length(ma) }
  }
  
  grad<- gradient(ar=ar,ma=ma,y=y)
  eps <- grad$eps
  der.eps <- grad$gradient
  
  J <- matrix(0, nrow = (p + q), ncol = (p + q))
  Upsilon <- matrix(0, nrow = p + q, ncol = n)
  sig2 <- mean(eps^2)
  J <- (2 / n) * der.eps %*% t(der.eps)
  if (kappa(J) < grand) matJ.inv <- solve(J) else {matJ.inv <- ginv(J)}
  
 
  
  Upsilon <- as.matrix(sapply(1:n, function(i) (2 * t(der.eps[, i]) * eps[i])))
  if (p+q==1){
    Upsilon<-t(Upsilon)
  }
  
  matI <- matXi(data = Upsilon, p = p, q = q)
  matOmega <- matJ.inv %*% matI %*% matJ.inv
  matvar.strong <- 2 * sig2 * matJ.inv
  ecart.type.Omega <- sqrt(diag(matOmega) / n)
  ecart.type.strong <- sqrt(diag(matvar.strong) / n)
  
  list(matJ = J, matI = matI, matJ.inv = matJ.inv, matOmega = matOmega, matvar.strong = matvar.strong, standard.dev.Omega = ecart.type.Omega, standard.dev.strong = ecart.type.strong, sig2 = sig2)
}
