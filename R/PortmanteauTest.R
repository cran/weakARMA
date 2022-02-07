#'Portmanteau tests for one lag.
#'
#' @description Computes Box-Pierce and Ljung-Box statistics for standard, modified and 
#' self-normalized test procedures.
#'
#' @param ar Vector of AR coefficients. If \code{NULL}, it is a MA process.
#' @param ma Vector of MA coefficients. If \code{NULL}, it is an AR process.
#' @param y Univariate time series.
#' @param h Integer for the chosen lag.
#' @param grad Gradient of the series from the function \link{gradient}. If \code{NULL} gradient will be computed.
#'
#' @importFrom matrixStats rowCumsums
#' @importFrom stats pchisq
#' @importFrom CompQuadForm imhof
#' @importFrom MASS ginv
#' 
#' @details Portmanteau statistics are generally used to  test the null hypothesis.
#'  H0 : \eqn{X_t} satisfies an ARMA(p,q) representation.
#'
#'  The Box-Pierce (BP) and Ljung-Box (LB) statistics, defined as follows, are
#'  based on the residual empirical autocorrelation. \deqn{Q_{m}^{BP} =
#'  n\sum_{h}^{m} \rho^{2}(h)} \deqn{Q_{m}^{LB} = n(n+2) \sum_{h}^{m}
#'  \frac{\rho^{2}(h)}{(n-h)}}
#'
#'  The standard test procedure consists in rejecting the null hypothesis of an
#'  ARMA(p,q) model if the statistic \eqn{Q_m > \chi^{2}(1-\alpha)} where
#'  \eqn{\chi^{2}(1-\alpha)} denotes the \eqn{(1-\alpha)}-quantile of a
#'  chi-squared distribution with m-(p+q) (where m > p + q) degrees of freedom. The
#'  two statistics have the same asymptotic distribution, but the LB statistic
#'  has the reputation of doing better for small or medium sized samples.
#'
#'  But the significance limits of the residual autocorrelation can be very
#'  different for an ARMA models with iid noise and ARMA models with only
#'  uncorrelated noise but dependant. The standard test is obtained under the
#'  stronger assumption that \eqn{\epsilon_{t}} is iid. So we give an another
#'  way to obtain the exact asymptotic distribution of the standard portmanteau
#'  statistics under the weak dependence assumptions.
#'
#'  Under H0, the statistics \eqn{Q_{m}^{BP}} and \eqn{Q_{m}^{LB}} converge in
#'  distribution as \eqn{n \rightarrow \infty}, to \deqn{Z_m(\xi_m) :=
#'  \sum_{i}^{m}\xi_{i,m} Z^{2}_i} where \eqn{\xi_m =
#'  (\xi_{1,m}',...,\xi_{m,m}')} is the eigenvalues vector of the asymptotic
#'  covariance matrix of the residual autocorrelations vector and
#'  \eqn{Z_{1},...,Z_{m}} are independent \eqn{\mathcal{N}(0,1)} variables.
#'  
#'
#'  So when the error process is a weak white noise, the asymptotic distribution
#'  \eqn{Q_{m}^{BP}} and \eqn{Q_{m}^{LB}} statistics is a weighted sum of
#'  chi-squared. The distribution of the quadratic form \eqn{Z_{m}(\xi_m)} can
#'  be computed using the algorithm by Imhof available here :
#'   \code{\link[CompQuadForm]{imhof}}
#'
#'  We propose an alternative method where we do not estimate an asymptotic
#'  covariance matrix. It is based on a self-normalization based approach to
#'  construct a new test-statistic which is asymptotically distribution-free
#'  under the null hypothesis.
#'
#'  The sample autocorrelation, at lag \code{h} take the form \eqn{\hat{\rho}(h) =
#'  \frac{\hat{\Gamma}(h)}{\hat{\Gamma}(0)}}. 
#'  Where \eqn{\hat{\Gamma}(h) = \frac{1}{n} \sum_{t=h+1}^n \hat{e}_t\hat{e}_{t-h}}.
#'  With \eqn{\hat{\Gamma}_m = (\hat{\Gamma}(1),...,\hat{\Gamma}(m))   }
#'  The vector of the first m sample autocorrelations is written \eqn{\hat{\rho}_m = (\hat{\rho}(1),...,\hat{\rho}(m))'}.
#'
#'  The normalization matrix is defined by \eqn{\hat{C}_{m} =
#'  \frac{1}{n^{2}}\sum_{t=1}^{n} \hat{S}_t \hat{S}_t'} where \eqn{\hat{S}_t = \sum_{j=1}^{t} (\hat{\Lambda} \hat{U}_{j} -
#'  \hat{\Gamma}_m)}.
#'
#'  The sample autocorrelations satisfy \eqn{Q_{m}^{SN}=n\hat{\sigma}^{4}\hat{\rho}_m '
#'  \hat{C}_m^{-1}\hat{\rho}_m \rightarrow U_{m}}.
#'  
#'  \eqn{\tilde{Q}_{m}^{SN} =
#'  n\hat{\sigma}^{4}\hat{\rho}_{m}' D_{n,m}^{1/2}\hat{C}_{m}^{-1} D_{n,m}^{1/2}\hat{\rho}_{m} \rightarrow U_{m} }
#'  reprensating respectively the version modified of Box-Pierce (BP) and
#'  Ljung-Box (LB) statistics. Where \eqn{D_{n,m} = \left(\begin{array}{ccc} \frac{n}{n-1} & & 0 \\
#'   & \ddots & \\
#'  0 & & \frac{n}{n-m}
#'  \end{array}\right)}.
#'  The critical values for \eqn{U_{m}} have been tabulated by Lobato.
#'  
#'
#'
#' @return A list including statistics and p-value:
#' 
#' \describe{
#'  \item{\code{Pm.BP}}{Standard portmanteau Box-Pierce statistics.}
#'  \item{\code{PvalBP}}{p-value corresponding at standard test where the
#'      asymptotic distribution is approximated by a chi-squared}
#'  \item{\code{PvalBP.Imhof}}{p-value corresponding at the exact asymptotic distribution
#'        of the standard portmanteau Box-Pierce statistics.}
#'  \item{\code{Pm.LB}}{Standard portmanteau Box-Pierce statistics. }
#'  \item{\code{PvalLB}}{p-value corresponding at standard test where the
#'      asymptotic distribution is approximated by a chi-squared. }
#'  \item{\code{PvalLB.Imhof}}{ p-value corresponding at the exact asymptotic distribution
#'        of the standard portmanteau Ljung-Box statistics.}
#'  \item{\code{LB.modSN }}{Ljung-Box statistic with the self-normalization method. }
#'  \item{\code{BP.modSN}}{Box-Pierce statistic with the self-normalization method.} }
#' 
#' @references Boubacar Maïnassara, Y. 2011, Multivariate portmanteau test for structural {VARMA} models
#' with uncorrelated but non-independent error terms \emph{Journal of Statistical Planning and Inference},
#' vol. 141, no. 8, pp. 2961-2975.
#' @references Boubacar Maïnassara, Y. and Saussereau, B. 2018, Diagnostic checking in multivariate {ARMA} models with 
#' dependent errors using normalized residual autocorrelations ,
#'  \emph{Journal of the American Statistical Association}, vol. 113, no. 524, pp. 1813-1827.
#' @references Francq, C., Roy, R. and Zakoïan, J.M. 2005, Diagnostic Checking in ARMA
#'  Models with Uncorrelated Errors, \emph{Journal of the American Statistical
#'  Association}, vol. 100, no. 470 pp. 532-544
#' @references Lobato, I.N. 2001, Testing that a dependant process is
#'  uncorrelated. J. Amer. Statist. Assos. 96, vol. 455, pp. 1066-1076.
#'
#' @seealso  \code{\link{portmanteauTest}} to obtain the statistics of all m
#'  lags.

portmanteauTest.h <- function(ar = NULL, ma = NULL, y, h, grad=NULL)
{
  grand<-1 / sqrt(.Machine$double.eps)  #nombre grand
  n <- length(y)  #nb d'obs
  h.max <- as.integer(min(10*log10(n), n-1))
  if (is.null(h)|| (h>h.max)) h <- as.integer(min(10*log10(n), n-1))
  h <- abs(floor(h))  #on veut un entier positif
  
  if (h>= n)return(NA) 
  
  if (is.null(ar))
  {
    if (is.null(ma)){
      p <- 0 ; q <- 0
    } else {
      p <- 0 ; q <- length(ma)
    }
  }
  else{
    if (is.null(ma)){
      q <- 0 ; p <- length(ar)}
    else{
      q <- length(ma) ; p <- length(ar)
    }
  } #valeurs p et q en fonction de ar et ma
  
  eps <- y  #Dans le cas ou p+q==0 eps=y, dans les autres eps sera remplace
  der.eps <- matrix(0, nrow = (p + q), ncol = n)
  upsilon <- matrix(0, nrow = (p + q), ncol = n)
  upsilon2 <- matrix(0, nrow = h, ncol = n)
  evec <- matrix(0, nrow = h, ncol = n)
  Prod <- array(0, dim = c(h, (p + q), n))
  J <- matrix(0, nrow = (p + q), ncol = (p + q))
  Phim <- array(0, dim = c(h, (p + q)))
  Upsilon <- array(0, dim = c(h + p + q, n))
  
  GAMMA <- matrix(0, nrow = h, ncol = (h + p + q))
  S_t <- matrix(0, nrow = h, ncol = n)
  matC_h <- matrix(0, nrow = h, ncol = h)
  Upsilon.reduit <- matrix(0, nrow = h, ncol = n)
  Upsilon_centre <- matrix(0, nrow = h, ncol = n)
  
  selec <- floor(n ^ ((1 / 3) - .Machine$double.eps))
  
  if ((p+q==0)){ #dans le cas ou p+q=0 on ne peut pas utiliser notre fonction gradient
    for (t in 2:n){
      for (j in 1:min(h,(t-1))){ #pour eviter iterations inutiles
        evec[j, t] <- eps[t - j] 
      }
    }
    upsilon2<-sapply(1:n, function(i) evec[,i]*eps[i]) #calcul ups2
    
  } else {
    if (is.null(grad)){
      grad <- gradient(ar = ar, ma = ma, y = y)
    }
    eps <- grad$eps
    der.eps <- grad$gradient
    
    upsilon<-sapply(1:n, function(i) eps[i]*der.eps[,i]) #calcul ups
    
    J <- 2 * der.eps %*% t(der.eps) / n #calcul J
    if (kappa(J) < grand) matJ.inv <- solve(J)
    else matJ.inv <- ginv(J) #calcul invJ
    
  }
  
  sig2 <- mean(eps^2)
  
  
  if (h==1){ #cas latence = 1
    if ((p+q==0)){ #cas bruit
      Upsilon <- sapply(1:n, function(i) as.vector(upsilon2[ i]))
      data <- as.vector(Upsilon) #estimation a besoin d un vect
      coef <- estimation(p = selec, y = data)
      phi <- 1 - sum(coef$ar)
      if (kappa(phi) < grand) phi.inv <- solve(phi)
      else phi.inv <- ginv(phi) #calcul inv phi
      
      Sigma_Gamma <- (phi.inv^2) * coef$sigma.carre #calcul sigma_gamma
    } else {
      if ((p+q==1)){ #cas AR ou MA
        for (t in 2:n){
          for (j in 1:min(h,(t-1))) {
            evec[j,t] <- eps[t - j]
          }
          Prod[ , ,t] <- evec[ ,t]%*%t(der.eps[ ,t])
        }  
        Phih <- mean(sapply(1:n, function(i) as.vector(Prod[ , ,i])))
        GAMMA <- c(1, Phih)
      } else { #cas ARMA ou AR,MA d ordre > 1
        for (t in 2:n){
          for (j in 1:min(h,(t-1))) {
            evec[j,t] <- eps[t - j]
          }
          Prod[ , ,t] <- evec[ ,t]%*%t(der.eps[ ,t])
        }  
        Phih <- matrix(rowMeans(sapply(1:n, function(i) as.vector(Prod[ , ,i]))), nrow = h, ncol = (p + q))
        GAMMA <- cbind(diag(1, h), Phih)
      }
    }
  } else {
    if ((p+q==0)){
      Upsilon<-sapply(1:n, function(i) as.vector(upsilon2[ ,i]))
      Sigma_Gamma <-matXi(Upsilon,p=h) #p=h pour ne pas que le if(p+q<=1) du matXi soit vrai et que  l'on obtienne une matrice [h x h]
    } else {
      for (t in 2:n){
        for (j in 1:min(h,(t-1))) {
          evec[j,t] <- eps[t - j]
        }
        Prod[ , ,t] <- evec[ ,t]%*%t(der.eps[ ,t])
      }  
      Phih <- matrix(rowMeans(sapply(1:n, function(i) as.vector(Prod[ , ,i]))), nrow = h, ncol = (p + q))
      GAMMA <- cbind(diag(1, h), Phih)
    }
  }
  
  upsilon2<-sapply(1:n, function(i) evec[,i]*eps[i])
  
  if ((p+q==0)){
    if (h==1){
      Upsilon.reduit <- sapply(1:n, function(i) as.vector((upsilon2[ i]))) #Upsilon.reduit pour harmoniser la sortie du if else 
    } else {
      Upsilon.reduit <- sapply(1:n, function(i) as.vector((upsilon2[ ,i]))) #Upsilon.reduit pour harmoniser la sortie du if else 
    }

  }else {
    if ((p+q==1)){
      if (h != 1) {
        Upsilon <- sapply(1:n, function(i)
          as.vector(c(upsilon2[ ,i], ((-2) * matJ.inv) %*% upsilon[ i])))
      } else {
        Upsilon <- sapply(1:n, function(i)
          as.vector(c(upsilon2[ i], ((-2) * matJ.inv) %*% upsilon[ i])))
      }
    } else {
      if (h != 1) {
        Upsilon <- sapply(1:n, function(i)
          as.vector(c(upsilon2[ ,i], ((-2) * matJ.inv) %*% upsilon[ ,i])))
      } else {
        Upsilon <- sapply(1:n, function(i)
          as.vector(c(upsilon2[ i], ((-2) * matJ.inv) %*% upsilon[ ,i])))
      }
    }
    
    
    Upsilon.reduit <- GAMMA %*% Upsilon
  }
  
  auto <- acf.gamma_m(ar = ar, ma = ma, y = y, h = h, e = eps)
  rho_h <- auto$rho_m
  gamma_h <- auto$gamma_m
  
  Upsilon.centre <- Upsilon.reduit-(matrix(t(gamma_h), nrow = h, ncol = n))
  
  S_t <- (rowCumsums(Upsilon.centre))
  
  for (i in 1:n){
    matC_h <- matC_h + (S_t[ ,i] %*% (t(S_t[ ,i]))) / (n^2)
  }
  
  if (kappa(matC_h) < grand)
    matC_h.inv <- solve(matC_h) else matC_h.inv <- ginv(matC_h)
  
 
  T_n <- matrix(0, nrow = h, ncol = h)
  for (i in 1:h){
    T_n[i, i] <- sqrt((n + 2)*(1 / (n-i)))
  }
  
  LB.modSN <- as.numeric( n * ((sig2)^2) * ( (t(rho_h) %*% T_n) %*% matC_h.inv %*% (T_n %*% rho_h))  )
  BP.modSN <- as.numeric( n * ((sig2)^2) * (t(rho_h) %*% matC_h.inv %*% (rho_h))   )
  
  if ((p+q!=0)){
    matI <- matXi(Upsilon,p=(p+h),q=q) #On calcule une matrice d'estimation des param pas la matI, p+m permet d'eviter if(p+q==1) dans le cas AR(1) ou MA(1)
                                        #car on veut une matrice dim (p+q+h x p+q+h)
    
    Sigma_gam <- matI[1:h, 1:h]
    Sigma_theta <- matI[(h + 1):(h + p + q), (h + 1):(h + p + q)]
    Sigma_gam.theta <- matI[1:h, (h + 1):(h + p + q)]
    Sigma_theta.gam <- matI[(h + 1):(h + p + q), 1:h]
    
    Sigma_Gamma <- Sigma_gam + Phih %*% Sigma_theta %*% t(Phih) + Phih %*% Sigma_theta.gam + Sigma_gam.theta %*% t(Phih)
    
  }
  
  Sigma_rho <- as.matrix(Sigma_Gamma / ((sig2)^2))
  
  Lambda <- eigen(Sigma_rho, symmetric = TRUE)$values
  
  Pm.LB <- as.numeric(n * (t(T_n %*% rho_h) %*% (T_n %*% rho_h))  )
  Pm.BP <- as.numeric(n * (t(rho_h) %*% ((rho_h)))   )
  
  if (h - (p + q) <= 0) {
    PvalLB <- NA
    PvalBP <- NA
  } else {
    PvalLB <- as.numeric( 1 - pchisq(Pm.LB, df = (h - (p + q)))  )
    PvalBP <- as.numeric( 1 - pchisq(Pm.BP, df = (h - (p + q)))  )
  }
  
  PvalLB.Imhof <- imhof(Pm.LB, Lambda)$Qq
  PvalBP.Imhof <- imhof(Pm.BP, Lambda)$Qq
  
  list(Omegah = Sigma_rho, Pm.BP = Pm.BP, PvalBP = PvalBP, PvalBP.Imhof = PvalBP.Imhof,
      Pm.LB = Pm.LB, PvalLB = PvalLB, PvalLB.Imhof = PvalLB.Imhof, LB.modSN = LB.modSN, BP.modSN = BP.modSN)
  
}

#'Portmanteau tests
#'
#' @description Realizes portmanteau tests of the first m lags, this function uses \code{\link{portmanteauTest.h}}
#' for h in 1:m.
#'
#' @param ar Vector of AR coefficients. If \code{NULL}, it is a MA process.
#' @param ma Vector of MA coefficients. If \code{NULL}, it is an AR process.
#' @param y Univariate time series.
#' @param m Integer for the lag.
#'
#' @return A list of vectors of length \code{m}, corresponding to statistics and p-value for each lag, 
#' for standard, modified and self-normalized Ljung-Box and Box-Pierce methods.
#'
#'  
#' @export
#'  
#' @examples 
#'  est<-estimation(p = 1, q = 1, y = CAC40return.sq)
#'  \donttest{portmanteauTest(ar = est$ar, ma = est$ma, y = CAC40return.sq, m = 20)}
#'
#' @references Boubacar Maïnassara, Y. 2011, Multivariate portmanteau test for structural {VARMA} models
#' with uncorrelated but non-independent error terms \emph{Journal of Statistical Planning and Inference},
#' vol. 141, no. 8, pp. 2961-2975.
#' @references Boubacar Maïnassara, Y. and Saussereau, B. 2018, Diagnostic checking in multivariate {ARMA} models with 
#' dependent errors using normalized residual autocorrelations ,
#'  \emph{Journal of the American Statistical Association}, vol. 113, no. 524, pp. 1813-1827.
#' @references Francq, C., Roy, R. and Zakoïan, J.M. 2005, Diagnostic Checking in ARMA
#'  Models with Uncorrelated Errors, \emph{Journal of the American Statistical
#'  Association}, vol. 100, no. 470, pp. 532-544.
#' 
#' 
#' @seealso \code{\link{portmanteauTest.h}} to obtain statistics for only one h lag.

portmanteauTest <- function(ar = NULL, ma = NULL, y, m = NULL)
{

  n <- length(y)
  m.max <- as.integer(min(10*log10(n), n-1))
  
  if (is.null(m)||m>m.max) m <- m.max
  
  Pm.BP <- rep(0, m)
  names(Pm.BP) <- paste("m", 1:m, sep = " = ")
  PvalBP <- rep(0, m)
  names(PvalBP) <- paste("m", 1:m, sep = " = ")
  PvalBP.Imhof <- rep(0, m)
  names(PvalBP.Imhof) <- paste("m", 1:m, sep = " = ")
  Pm.LB <- rep(0, m)
  names(Pm.LB) <- paste("m", 1:m, sep = " = ")
  PvalLB <- rep(0, m)
  names(PvalLB) <- paste("m", 1:m, sep = " = ")
  PvalLB.Imhof <- rep(0, m)
  names(PvalLB.Imhof) <- paste("m", 1:m, sep = " = ")
  LB.modSN <- rep(0, m)
  names(LB.modSN) <- paste("m", 1:m, sep = " = ")
  BP.modSN <- rep(0, m)
  names(BP.modSN) <- paste("m", 1:m, sep = " = ")
  
  if (is.null(ar) & is.null(ma)){
    grad<-NULL
  } else {
    grad<-gradient(ar = ar, ma = ma, y = y)
  }
  
  
  for(i in 1:m){
    res <- portmanteauTest.h(ar = ar, ma = ma, y = y, h = i, grad = grad)
    
    Pm.BP[i] <- round( res$Pm.BP, digits = 6)
    PvalBP[i] <- round(res$PvalBP, digits = 6)
    PvalBP.Imhof[i] <- round(res$PvalBP.Imhof, digits = 6)
    Pm.LB[i] <- round(res$Pm.LB, digits = 6)
    PvalLB[i] <- round(res$PvalLB, digits = 6)
    PvalLB.Imhof[i] <- round(res$PvalLB.Imhof, digits = 6)
    LB.modSN[i] <- round(res$LB.modSN, digits = 6)
    BP.modSN[i] <- round(res$BP.modSN, digits = 6)
    
  }
  
  list(Pm.BP = Pm.BP, PvalBP = PvalBP, PvalBP.Imhof = PvalBP.Imhof,
       Pm.LB = Pm.LB, PvalLB = PvalLB, PvalLB.Imhof = PvalLB.Imhof, LB.modSN = LB.modSN, BP.modSN = BP.modSN)
  
}

resultat.h <- function (ar = NULL, ma = NULL, y, m = NULL, eps)
{
  
  grand<-1 / sqrt(.Machine$double.eps)
  n <- length(y)
  m <- abs(m)
  selec <- min(floor(n ^ ((1 / 3) - .Machine$double.eps)),5)
  if (m >= n) return(NA)
  
 
  
  auto <- acf.univ(ar = ar, ma = ma, y = y, h = m, e = eps)
  rho_m <- auto$autocor
  gamma_m <- auto$autocov
  
  ee.h <- matrix(0, nrow = 1, ncol = n)
  upsilon.h <- matrix(0, nrow = 1, ncol = n)
  Upsilon.reduit.h <- matrix(0, nrow = 1, ncol = n)
  Upsilon_centre.h <- matrix(0, nrow = 1, ncol = n)
  
  
  if (is.null(ar) & is.null(ma))
  {
    eps <- y
    sig2 <- mean(eps^2)
    
    
    for (t in (m+1):n){
      upsilon.h[ ,t] <- eps[t - m] * eps[t]
    }
    
    Upsilon.h <- sapply(1:n, function(i) as.vector(upsilon.h[1,i]))
    Upsilon.reduit.h <- Upsilon.h
    Upsilon_centre.h <- as.vector(Upsilon.reduit.h-(array(t(gamma_m), dim = c(1, n))))
    S_t.h <- cumsum(Upsilon_centre.h) #vecteur de la somme cumulée 
    matC_m.h <- 0
    for (i in 1:n) matC_m.h <- matC_m.h + (S_t.h[i] * (t(S_t.h[i]))) / (n ^ 2)
    matC_m.h <- matC_m.h / (sig2^2)
    if (kappa(matC_m.h) < grand) matC_m.inv.h <- solve(matC_m.h)
    else matC_m.inv.h <- matC_m.h
    
    data.h <- Upsilon.h[1:n]
 
    coef <- estimation(p = selec, y = data.h)
    phi <- 1 - sum(coef$ar)
    if (kappa(phi) < grand) phi.inv <- solve(phi) else phi.inv <- phi
    matI.h <- as.numeric((phi.inv^2)) * coef$sigma.carre
    Sigma_rho.h <- as.numeric(matI.h / ((sig2)^2))
    
    return(list(Sigma_rho.h = Sigma_rho.h, eps = eps, matC_m.inv.h = matC_m.inv.h, matC_m.h = matC_m.h))
  } else {
    
      if (is.null(ma))
      { q <- 0 ; p <- length(ar)}
      else {
        if (is.null(ar))
        { p <- 0 ; q <- length(ma)}
        else
        {q <- length(ma) ; p <- length(ar)}
      }
  
    
    der.eps <- matrix(0, nrow = (p + q), ncol = n)
    upsilon <- matrix(0, nrow = (p + q), ncol = n)
    Prod.h <- array(0, dim = c(1, (p + q), n))
    J <- matrix(0, nrow = (p + q), ncol = (p + q))
    Phim.h <- array(0, dim = c(1, (p + q)))
    Upsilon.h <- matrix(0, nrow = (1 + p + q ), ncol = n)
    GAMMA.h <- matrix(0, nrow = 1, ncol = (p + q + 1))
    S_t.h <- array(0, dim = c(1, n))
    
    
    
    
    grad <- gradient(ar = ar, ma = ma, y = y)
    eps <- grad$eps
    der.eps <- grad$gradient
    sig2 <- mean(eps^2)
    
    J <- 2*der.eps %*% t(der.eps) / n
    if (kappa(J)< grand )matJ.inv <- solve(J) else {matJ.inv <- ginv(J)}
    
  
    
    
    upsilon <- sapply(1:n, function (i) eps[i]*der.eps[ ,i]) 
    for (t in (m+1):n){
      upsilon.h[t] <- eps[t - m]*eps[t]
      if ((p+q==1))
        Prod.h[1, ,t] <- eps[t-m] * t(der.eps[ ,t])
      else  Prod.h[ , ,t] <- eps[t-m] %*% t(der.eps[ ,t])
    }
    

    if ((p+q)==1)
    {
      
      Phim.h <- mean(sapply(1:n, function(i) Prod.h[ , ,i])) 
      Upsilon.h <- sapply(1:n, function(i) c(upsilon.h[1,i], ((- 2) * matJ.inv) * upsilon[i]))
      GAMMA.h <- cbind(1, Phim.h)
    }
    else{
      Phim.h <- matrix(rowMeans(sapply(1:n, function(i) Prod.h[ , ,i])), nrow = 1, ncol = (p + q))
      Upsilon.h <- as.matrix( sapply(1:n, function(i) c(upsilon.h[1,i], ((-2) * matJ.inv) %*% upsilon[ ,i])))
      GAMMA.h <- cbind(1, Phim.h)
    }
    
    Upsilon.reduit.h <- GAMMA.h %*% Upsilon.h
    Upsilon_centre.h <- as.vector(Upsilon.reduit.h - (array(t(gamma_m), dim = c(1, n))))
    S_t.h <- cumsum(Upsilon_centre.h)
    matC_m.h <- 0
    for (i in 1:n) matC_m.h <- matC_m.h + (S_t.h[i]*(t(S_t.h[i]))) / (n^2)
    matC_m.h <- matC_m.h / (sig2^2)
    if (kappa(matC_m.h)< grand) matC_m.inv.h <- solve(matC_m.h)
    else {matC_m.inv.h <- ginv(matC_m.h)}
    
    matXI.h <- matXi(Upsilon.h,(p+1),q) 
    
    Sigma_gam.h <- matXI.h[1, 1]
    Sigma_theta.h <- matXI.h[2:(p + q + 1), 2:(p + q + 1)]
    Sigma_gam.theta.h <- matXI.h[1, 2:(p + q + 1)]
    Sigma_theta.gam.h <- matXI.h[2:(p + q + 1), 1]
    
    Sigma_Gamma.h <- Sigma_gam.h + Phim.h %*% Sigma_theta.h %*% t(Phim.h) + Phim.h %*% Sigma_theta.gam.h +
      Sigma_gam.theta.h %*% t(Phim.h)
    Sigma_rho.h <- as.numeric(Sigma_Gamma.h / ((sig2)^2))
    
    return(list(Sigma_rho.h = Sigma_rho.h, eps = eps,
                matC_m.inv.h = matC_m.inv.h, matC_m.h = matC_m.h))
    
  }
  
}
