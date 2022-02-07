residual <- function(ar = NULL, ma = NULL, y, centre=TRUE)
{
  n <- length(y)
  e <- y
  mu<-0
  if (centre==FALSE){
    mu<-mean(y)
  }
  ybis<-y-mu
  if (is.null(ar) && (is.null(ma)))
    e <- y
  else{
    if (is.null(ma))
      for (t in 2:n){ e[t] <- ybis[t]- sum(ar[1:min((t-1),length(ar))]*ybis[(t-1):max(1,(t - length(ar)))]) }
    else
      if ((is.null(ar)))
        for (t in 2:n)
          e[t] <- ybis[t]+ sum(ma[1:min((t-1),length(ma))]*e[(t-1):max(1,(t -length(ma)))])
      else{
        for (t in 2:n) { 
          e[t] <- ybis[t]- sum(ar[1:min((t-1),length(ar))]*ybis[(t-1):max(1,(t - length(ar)))]) + sum(ma[1:min((t-1),length(ma))]*e[(t-1):max(1,(t - length(ma)))])
        }
      }
  }
  return(e)
}
