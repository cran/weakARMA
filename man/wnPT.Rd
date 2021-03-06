% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WhiteNoise.R
\name{wnPT}
\alias{wnPT}
\title{Weak white noise}
\usage{
wnPT(n, sigma = 1, k = 1, ninit = 100)
}
\arguments{
\item{n}{Number of observations.}

\item{sigma}{Standard deviation.}

\item{k}{Integer corresponding to the number of past observation will be used.}

\item{ninit}{Length of 'burn-in' period.}
}
\value{
Vector of size \code{n} containing  a nonlinear sequence \eqn{X_{i}} such as
  \eqn{X_{i} = Z_{i}Z_{i-1}...Z_{i-k}} , where \eqn{Z_{i}} is a sequence of iid
    random variables mean-zero random variable with variance \eqn{\sigma^2}.
}
\description{
Simulates an uncorrelated but dependant noise process.
}
\examples{
wnPT(100)
wnPT(100, sigma = 1, k = 1)
wnPT(100, k = 0) #strong noise
}
\references{
Romano, J. and Thombs, L. 1996, Inference for autocorrelation under weak assumptions,
 \emph{Journal of the American Statistical Association}, vol. 91, no. 434, pp. 590-600
}
\seealso{
\code{\link{wnRT}}, \code{\link{wnPT_SQ}}, \code{\link{simGARCH}}
}
