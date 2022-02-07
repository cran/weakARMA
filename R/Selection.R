#'Selection of ARMA models
#'
#' @description Identifies the orders p and q of an ARMA model according to several
#'  information criteria.
#'
#' @param data  Univariate time series.
#' @param P Integer for the maximum lag order of autoregressive component.
#' @param Q Integer for the maximum lag order of moving-average component.
#' @param c Real number >1 needed to compute Hannan-Quinn information criterion.
#'
#'
#'
#' @return  A list of the different criteria, each item contains the matrix of the 
#'  computed value for the different model and the selected order with this criterion 
#'  (corresponding to the minimum value in the previous matrix).
#'
#' @details The fitted model which is favored is the one corresponding to the
#'  minimum value of the criterion. The most popular criterion is the Akaike information
#'  criterion (\code{AIC}). This was designed to be an approximately unbiased
#'  estimator of a fitted model. For small sample or when  the number of fitted
#'  parameters is large, it is more appropriate to  manipulate a corrected AIC
#'  version (\code{AICc}) which is more nearly unbiased. But these two criteria
#'  are inconsistent for model orders selection. If you want to use a consistent
#'  criterion, it is possible to take the Bayesian information criterion
#'  (\code{BIC}) or the Hannan-Quinn information  criteria (\code{HQ}).
#'
#'  For the weak ARMA, i.e under the assumption that the errors are uncorrelated
#'  but not necessarily independant, modified criteria has been adapted :
#'  \code{AICm}, \code{AICcm}, \code{BICm}, \code{HQm}.
#'
#'  The criteria definitions are the following :
#'
#'  \deqn{AIC = n\log(\sigma^{2}) + 2(p + q)}
#'
#'  \deqn{AICm = n\log(\sigma^{2}) + \frac{Tr(IJ^{-1})}{\sigma^2}}
#'
#'  \deqn{AICc = n\log(\sigma^{2}) + n + \frac{n}{(n-(p + q + 1))} 2(p + q)}
#'
#'  \deqn{AICcm = n\log(\sigma^{2}) + \frac{n^{2}}{(n-(p + q + 1))}  + \frac{n}{(2(n-(p + q + 1)))} \frac{Tr(IJ^{-1})}{\sigma^2}}
#'
#'  \deqn{BIC = n\log(\sigma^{2}) + (p + q)log(n)}
#'
#'  \deqn{BICm = n\log(\sigma^{2}) + \frac{1}{2} \frac{Tr(IJ^{-1})}{\sigma^2}log(n)}
#'
#'  \deqn{HQ = n\log(\sigma^{2}) + 2c(p + q)log(log(n))}
#'
#'  \deqn{HQm = n\log(\sigma^{2}) + c\frac{Tr(IJ^{-1})}{\sigma^2}log(log(n))}
#'
#'
#' @export 
#'
#' @examples \donttest{ARMA.selec (CAC40return.sq, P = 3, Q = 3)}
#'
#' @references Boubacar Maïnassara, Y. 2012, Selection of weak VARMA models by
#'  modified Akaike's information  criteria, \emph{Journal of Time Series
#'  Analysis}, vol. 33, no. 1, pp. 121-130
#' @references Boubacar Maïnassara, Y. and Kokonendji, C. C. 2016, Modified Schwarz
#'  and Hannan-Quin information criteria for weak VARMA models, \emph{Stat
#'  Inference Stoch Process}, vol. 19, no. 2, pp. 199-217

ARMA.selec <- function(data, P, Q, c = 2)
{
  I <- array(0, c(P + Q + 1, P + Q + 1, (P+1)*(Q+1)))
  J.inv <- array(0, c(P + Q + 1, P + Q + 1, (P+1)*(Q+1)))
  sigma <- rep(0, (P+1)*(Q+1))
  
  i <- 1
  for(p in 0:P){
    for(q in 0:Q){
      if ((p+q)!=0){
        para.estim <- estimation(p = p, q = q, y = data)
        if (p!=0){
          estimateur.ar <- as.vector(para.estim$ar)
        } else {
          estimateur.ar <- NULL
        }
        if (q!=0){
          estimateur.ma <- as.vector(para.estim$ma)
        } else {
          estimateur.ma <- NULL
        }
        omeg <- omega(ar = estimateur.ar, ma = estimateur.ma, y = data)
        I[1:(p + q), 1:(p + q), i] <- as.matrix(omeg$matI)
        J.inv[1:(p + q), 1:(p + q), i] <- omeg$matJ.inv
        sigma[i] <- omeg$sig2
      } 
      i <- i + 1
    }
  }
  
  AIC <- critere.ARMA(data = data, P = P, Q = Q, func = AICb, c = c, I = I, J.inv = J.inv, sigma = sigma, critere.bruit = 2)
  AICc <- critere.ARMA(data = data, P = P, Q = Q, func = AICc, c = c, I = I, J.inv = J.inv, sigma = sigma, critere.bruit = 1)
  AICcm <- critere.ARMA(data = data, P = P, Q = Q, func = AICcm, c = c, I = I, J.inv = J.inv, sigma = sigma, critere.bruit = 1)
  AICm <- critere.ARMA(data = data, P = P, Q = Q, func = AICm, c = c, I = I, J.inv = J.inv, sigma = sigma, critere.bruit = 2)
  BIC <- critere.ARMA(data = data, P = P, Q = Q, func = BICb, c = c, I = I, J.inv = J.inv, sigma = sigma, critere.bruit = 2)
  BICm <- critere.ARMA(data = data, P = P, Q = Q, func = BICm, c = c, I = I, J.inv = J.inv, sigma = sigma, critere.bruit = 2)
  HQ <- critere.ARMA(data = data, P = P, Q = Q, func = HQ, c = c, I = I, J.inv = J.inv, sigma = sigma, critere.bruit = 2)
  HQm <- critere.ARMA(data = data, P = P, Q = Q, func = HQm, c = c, I = I, J.inv = J.inv, sigma = sigma, critere.bruit = 2)
  
  res <- list(AIC = AIC, AICm = AICm, AICc = AICc, AICcm = AICcm, BIC = BIC, BICm = BICm, HQ = HQ, HQm = HQm)
}

critere.ARMA <- function(data, P, Q, func, c = 2, I, J.inv, sigma, critere.bruit)
{
  n <- length(data)
  critere <- matrix(0, nrow = (P+1), ncol = (Q+1))
  
  
  
  i <- 1
  for(p in 0:P){
    for(q in 0:Q){
      if ((p+q)==0){
        switch (critere.bruit,
                critere[1,1] <- n*log(mean(data^2)) + n,
                critere[1,1] <- n*log(mean(data^2))
        )
        
      }else {
        critere[(p+1), (q+1)] <- func(n = n, p = p, q = q, sigma = sigma[i], I = I[1:(p + q), 1:(p + q), i], J.inv = J.inv[1:(p + q), 1:(p + q), i], c = c)
      }
      i <- i + 1
    }
  }
  
  rownames(critere)<-(seq(1:(P+1))-1)
  colnames(critere)<-(seq(1:(Q+1))-1)
  
  min.ind <- which(critere == min(critere), arr.ind = T)
  p <- (min.ind[, 1]-1) ; q <- (min.ind[, 2]-1)
  
  
  
  return(list(critere = critere, p = p, q = q))
}

