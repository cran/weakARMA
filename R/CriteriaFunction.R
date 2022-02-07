AICb <- function(n, p = NULL, q = NULL, sigma2, I, J.inv, c = 2)
{
  if (is.null(p)) p = 0
  if (is.null(q)) q = 0
  return(n*log(sigma2) + 2*(p + q))
}

AICc <- function(n, p = NULL, q = NULL, sigma2, I, J.inv, c = 2)
{
  if (is.null(p)) p = 0
  if (is.null(q)) q = 0
  return(n*log(sigma2) + n + (n / (n-(p + q))) * 2*(p + q))
}

AICcm <- function(n, p = NULL, q = NULL, sigma2, I, J.inv, c = 2)
{
  if (is.null(p)) p = 0
  if (is.null(q)) q = 0
  return(n*log(sigma2) + (n^2 / (n-(p + q ))) + (n / (2*(n-(p + q ))))*sum(diag((I%*%J.inv))/sigma2))
}

AICm <- function(n, p = NULL, q = NULL, sigma2, I, J.inv, c = 2)
{
  if (is.null(p)) p = 0
  if (is.null(q)) q = 0
  return(n*log(sigma2) + sum(diag((I%*%J.inv))/sigma2))
}

BICb <- function(n, p = NULL, q = NULL, sigma2, I, J.inv, c = 2)
{
  if (is.null(p)) p = 0
  if (is.null(q)) q = 0
  return(n*log(sigma2) + (p + q)*log(n))
}

BICm <- function(n, p = NULL, q = NULL, sigma2, I, J.inv, c = 2)
{
  if (is.null(p)) p = 0
  if (is.null(q)) q = 0
  return(n*log(sigma2) + 0.5*sum(diag((I%*%J.inv))/sigma2)*log(n))
}

HQ <- function(n, p = NULL, q = NULL, sigma2, I, J.inv, c = 2)
{
  if (is.null(p)) p = 0
  if (is.null(q)) q = 0
  return(n*log(sigma2) + c*2*(p + q)*log(log(n)))
}

HQm <- function(n, p = NULL, q = NULL, sigma2, I, J.inv, c = 2)
{
  return(n*log(sigma2) + c*sum(diag((I%*%J.inv))/sigma2)*log(log(n)))
}
