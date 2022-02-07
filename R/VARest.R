#' Estimation of VAR(p) model
#'
#' @param x Matrix of dimension (n,p+q).
#' @param p Integer for the lag order.
#'
#'
#' @return A list containing:
#' \describe{
#'    \item{\code{ac}}{Coefficients data matrix.}
#'    \item{\code{p}}{Integer of the lag order.}
#'    \item{\code{k}}{Dimension of the VAR.}
#'    \item{\code{res}}{Matrix of residuals.}
#' }
#' 
#' @description Estimates the coefficients of a VAR(p) model. Used in \code{\link{matXi}}.
#' 
#' 
#' @export
#'

VARest <- function (x, p) 
{
  n <- nrow(x)
  k <- ncol(x)
  y <- t(x[(p + 1):n, ])
  z0 <- t(x[n:1, ])
  z <- matrix(1, nrow = k * p, ncol = 1)
  ac <- array(1,dim = c(k,k,p))
  for (i in n:(p + 1)) {
    m <- t(t(as.vector(matrix(z0[, (i - p + 1):i]))))
    z <- cbind(z, m)
  }
 
  z <- z[, 2:(n - p + 1), drop = FALSE]

  if (kappa(z)<(1 / sqrt(.Machine$double.eps))){
    b <- tcrossprod(y, z) %*% solve(tcrossprod(z))
  }else {
    b <- tcrossprod(y, z) %*% ginv(tcrossprod(z))
  }
  e <- y - b %*% z

  for (i in 1:p){
    ac[,,i]<-b[,(1+(k*(i-1))):(k*i)]
  }
  
  return(list(ac=ac, p=p, k=k, res=t(e)))
}
