#' Estimation of Fisher information matrix I
#'  
#' @description Uses a consistent estimator of the matrix I based on an autoregressive spectral estimator.
#' 
#' @param data Matrix of dimension (p+q,n).
#' @param p Dimension of AR estimate coefficients.
#' @param q Dimension of MA estimate coefficients.
#' 
#' @importFrom stats resid
#'
#' @return Estimate Fisher information matrix \eqn{I = 
#' \sum_{h=-\infty}^{+\infty} cov(2e_t \nabla e_t, 2e_{t-h} \nabla e_{t-h})} where \eqn{\nabla e_t} 
#' denotes the gradient of the residuals.
#' 
#' @export
#'
#' @references Berk, Kenneth N. 1974, Consistent autoregressive spectral estimates,
#'  \emph{The Annals of Statistics}, vol. 2, pp. 489-502.
#' @references Boubacar Maïnassara, Y. and Francq, C. 2011, Estimating structural VARMA models with uncorrelated but 
#' non-independent error terms, \emph{Journal of Multivariate Analysis}, vol. 102, no. 3, pp. 496-505.  
#' @references Boubacar Mainassara, Y. and Carbon, M. and Francq, C. 2012, Computing and estimating information matrices 
#' of weak ARMA models \emph{Computational Statistics & Data Analysis}, vol. 56, no. 2, pp. 345-361.
#' 
#'  
#'
matXi <- function (data,p=0,q=0)
{
  if(length(data[,1])!=(p+q)){
    stop("Matrix data has not the right dimensions. Please use this code to obtain the right matrix:
        grad<- gradient(ar=ar,ma=ma,y=y)
        eps <- grad$eps
        der.eps <- grad$gradient
        data<-as.matrix(sapply(1:n, function(i) as.vector(2 * t(der.eps[, i]) * eps[i])))
        if ((length(ar)+length(ma))==1){
          data<-t(data)
        }
         with y your time serie and the correct ar and ma arguments.")
  }
  grand = 1 / sqrt(.Machine$double.eps)
  if ((p+q)<=1) {
    datab <- as.vector(data)
    n <- length(datab)
    selec <- floor(n^((1 / 3) - .Machine$double.eps))
    selec<-min(selec,5)
    coef <- estimation(p = selec, y = datab)
    phi <- 1 - sum(coef$ar)
    if (kappa(phi) < grand) phi.inv <- solve(phi) else phi.inv <- phi
    matXi <- (phi.inv^2) * coef$sigma.carre
  }else {
    data <- as.data.frame(t(data))
    n <- length(data[,1])
    selec <- floor(n^((1 / 3) - .Machine$double.eps))
    selec<-min(selec,5)
    mod<- VARest(x=data,p=selec)
    p1 <- mod$p
    d <- mod$k
    phi <- diag(1, d)
    Ac<-mod$ac
    for (i in 1:p1)
      phi <- phi - Ac[,,i]
    
    if (kappa(phi) < grand) phi.inv <- solve(phi) else phi.inv <- ginv(phi)
    res <- mod$res
    sigmaU <- (t(res) %*% res) / (n-p1)
    matXi <- phi.inv %*% sigmaU %*% t(phi.inv)
    
  }
  return(matXi)
  
}
