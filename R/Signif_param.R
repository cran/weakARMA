#' Computes the parameters significance
#' 
#' @description Computes a matrix with estimated coefficient and their significance.
#'
#' @param ar Vector of AR coefficients, if \code{NULL}, MA process.
#' @param ma Vector of MA coefficients, if \code{NULL}, AR process.
#' @param p Order of AR, if \code{NULL} MA process.
#' @param q Order of MA, if \code{NULL} AR process.
#' @param y Univariate time series.
#' @param sd.strong Standard error of time series in the strong case computed in \code{\link[weakARMA]{omega}}, if not provided the function will compute it.  
#' @param sd.weak Standard error of time series in the weak case computed in \code{\link[weakARMA]{omega}}, if not provided the function will compute it.
#' @param meanparam If \eqn{\mu} of the time series needs to be computed.
#' @param mu Value of \eqn{\mu}, if it is known and if the \code{meanparam} is \code{TRUE}. If not known the function will compute it. 
#'
#' @return Matrix of the estimate coefficient with their significance.
#' \describe{
#'     \item{\code{coef}}{Estimation of each coefficient.}
#'     \item{\code{sd}}{Standard deviation in each case.}
#'     \item{\code{t-ratio}}{T-ratio corresponding to each coefficient.}
#'     \item{\code{signif}}{Significance of each parameter. 
#'     Must be small, if not the parameter is not significant.}
#' }
#' 
#' @importFrom stats pt sd
#' 
#' @export
#'
#' @details The function needs at least one pair between: ar and/or ma, or p and/or q to be executed. It will be faster with all the parameters provided.
#' 
#' @examples 
#' \donttest{signifparam(p = 1, q = 2, y = CAC40return.sq)} #The last parameter is not significant.
#' \donttest{signifparam(p = 1, q = 1, y = CAC40return.sq)} #All the parameters are significant.
#' 

signifparam<-function(ar = NULL, ma = NULL, p = NULL, q = NULL, y, sd.strong = NULL, sd.weak = NULL, meanparam = TRUE, mu = NULL){
  if (is.null(p) & is.null(q) & is.null(ar) & is.null(ma)){
    stop("Impossible to compute without vectors ar and/or ma, or without order p and/or q")
  }
  if (is.null(p)){
    p<-length(ar)
  }
  if (is.null(q)){
    q<-length(ma)
  }
  
  n<-length(y)
  pValWeak<-numeric(p+q)
  pValStr<-numeric(p+q)
  names<-vector(length = (p+q))
  if ((is.null(ar) & is.null(ma))||(meanparam==TRUE & is.null(mu))){
    est<-estimation(p = p, q = q, y = y, meanparam = meanparam)
    mu<-est$mu
    ar<-est$ar
    ma<-est$ma
    sigma<-est$sigma.carre
  }
  if (is.null(sd.strong) || is.null(sd.weak)){
    om<-omega(ar = ar, ma = ma, y = y)
    sd.strong <- om$standard.dev.strong
    sd.weak <- om$standard.dev.Omega
  }
  tWeak<-c(ar,ma)/sd.weak
  tStr<-c(ar,ma)/sd.strong
  for (i in 1:(p+q)){
      pValWeak[i]<-round(2*(1-pt(abs(tWeak[i]),df=n-(p+q))),4)
      pValStr[i]<-round(2*(1-pt(abs(tStr[i]),df=n-(p+q))),4)
  }
  if (meanparam == TRUE){
    tMuW<-sqrt(n)*(mu/sqrt((((1-sum(ma))/(1-sum(ar)))^2)*(sigma)))
    pValMuW<-round(2*(1-pt(abs(tMuW),df=(n))),4)
    

    strg<-matrix(data = c(mu,t(c(ar,ma)),sqrt((((1-sum(ma))/(1-sum(ar)))^2)*(sigma/n)),t(sd.strong),mu/sqrt((((1-sum(ma))/(1-sum(ar)))^2)*(sigma/n)),t(tStr),pValMuW,t(pValStr)),nrow = (p+q+1),ncol=4)
    weak<-matrix(data = c(mu,t(c(ar,ma)),sqrt((((1-sum(ma))/(1-sum(ar)))^2)*(sigma/n)),t(sd.weak),mu/sqrt((((1-sum(ma))/(1-sum(ar)))^2)*(sigma/n)),t(tWeak),pValMuW,t(pValWeak)),nrow = (p+q+1),ncol=4)
    names[1]<-"mu"
    for ( i in 1:p){
      names[i+1]<-paste("alpha", i, sep=" ")
    }
    for ( i in 1:q){
      names[i+p+1]<-paste("beta", i, sep=" ")
    }
  }else {
    strg<-matrix(data = c(t(c(ar,ma)),t(tStr),t(sd.strong),t(pValStr)),nrow = (p+q),ncol=4)
    weak<-matrix(data = c(t(c(ar,ma)),t(tWeak),t(sd.weak),t(pValWeak)),nrow = (p+q),ncol=4)
    
    for ( i in 1:p){
      names[i]<-paste("alpha", i, sep=" ")
    }
    for ( i in 1:q){
      names[i+p]<-paste("beta", i, sep=" ")
    }
  }
  
  colnames(strg)<-c("coef","sd","t-ratio","signif")
  colnames(weak)<-c("coef","sd","t-ratio","signif")
  
  rownames(strg)<-names
  rownames(weak)<-names
  ret<-list(weak=weak,strong=strg)
  
  return(ret)
}

