#' @title Paris stock exchange return
#' @docType data
#' @description
#' This data set considers CAC40 return at the closure of the market from
#' March 2, 1990 to June 14, 2021.
#' 
#' @format A numerical vector with 7935 observations.
#' 
#' We computed every value from the dataset \code{\link{CAC40}}
#' with the following code:
#' 
#' \preformatted{
#'  cac<-CAC40;
#'  n<-length(cac);
#'  rend<-rep(0,n);
#'  rend[2:n]<-(log(cac[2:n]/cac[1:(n-1)])*100);
#'  CAC40return<-rend[2:n]
#'  }
#'  
#' @seealso \code{\link{CAC40}} and \code{\link{CAC40return.sq}}
#'  
"CAC40return"