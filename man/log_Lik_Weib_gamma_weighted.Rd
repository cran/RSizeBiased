\name{log_Lik_Weib_gamma_weighted}
\alias{log_Lik_Weib_gamma_weighted}


\title{Log likelihood function for the weighted gamma or Weibull distributions.}

\description{Calculates the log-likelihood function of the weighted gamma or Weibull (depends on user input)  distribution.}

\usage{log_Lik_Weib_gamma_weighted(TRpar,datain,r,dist)}

\arguments{
  \item{TRpar }{A vector of length 2, containing the shape and scale parameters of the distribution.}
  \item{datain }{The available sample points.}
   \item{r }{The size (order) of the distribution. The special cases \eqn{r=1,2,3} correspond to  length, area, volume biased samples respectively and are the most frequently encountered in practice. The case \eqn{r=0} corresponds to random samples from the Gamma distribution.}
    \item{dist }{Character switch, enables the choice of distribution: type "weib" for the Weibull or "gamma" for the gamma distribution.}
}

\details{The  log likelihood function of the weighted gamma  distribution is defined  by

\deqn{ \log L = \sum_{i=1}^n log f_r(X_i; \theta) }

where \eqn{f_r(x; \theta)} is the density of the \eqn{r-}size biased gamma distribution. Setting \eqn{r=0} corresponds to the log likelihood of the Gamma distribution.

In the case of Weibull, the  log likelihood   is defined  by

\deqn{ \log L = \sum_{i=1}^n log f_r(X_i; \theta) }

where \eqn{f_r(x; \theta)} is the density of the \eqn{r-}size biased Weibull distribution. Setting \eqn{r=0} corresponds to the log likelihood of the Weibull distribution.

}

\value{A scalar, the result of the log likelihood calculation.}

\references{ Economou et. al. (2021). Hypothesis testing for the population mean and variance based on r-size biased samples, under review.}

\author{
Polychronis Economou

R implementation and documentation: Polychronis Economou  <peconom@upatras.gr>
  }


%\seealso{ \code{\link{log_Lik_gam}} }

\examples{

#Log-likelihood for the gamma distribution for true parms=(2,3), r=0:
log_Lik_Weib_gamma_weighted(c(2,3), rgamma(100, shape=2, scale=3), 0, "gamma")
#Log-likelihood for the Weibull distribution for true parms=(2,3), r=0:
log_Lik_Weib_gamma_weighted(c(2,3), rweibull(100, shape=2, scale=3), 0, "weib")
}
