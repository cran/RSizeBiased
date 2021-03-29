\name{s11.s22}
\alias{s11.s22}


\title{Variance estimates for test statistics \eqn{\zeta_{n,r}^i, i=1,2} specifically for the Weibull and gamma distributions.}

\description{Variance estimates for test statistics \eqn{\zeta_{n,r}^i, i=1,2} specifically for the Weibull and gamma distributions.}

\usage{s11.s22(TRpar,r,sgg,dist)}


\arguments{
\item{TRpar }{A vector of length 2, containing the shape and scale parameters of the Weibull distribution.}
\item{r }{The size (order) of the distribution. The special cases \eqn{r=1,2,3} correspond to  length, area, volume biased samples respectively and are the most frequently encountered in practice. The case \eqn{r=0} corresponds to random samples from the underlying distribution.}
\item{sgg }{ Character switch ("s11" or "s22"), enables choosing between the s11 and s22 options}
\item{dist }{Character switch, enables the choice of distribution: type "weib" for the Weibull or "gamma" for the gamma distribution.}
}



\details{Provided that \eqn{\mu_r, r=1, 2, \dots } is the \eqn{r}th moment of the Weibull or the Gamma distribution, then

\deqn{ \sigma_{1,r}^2 = \mu_r (\mu_{2-r}) - 2 \mu_1 \mu_{1-r} + \mu_1^2 \mu_{-r}}

and

\deqn{ \sigma_{2,r}^2 = -4\mu_r \bigl ( 2\mu_{1}^2 - \mu_2) - 2) \mu_1 \mu_{1-r} +  (2\mu_1^2 - \mu_{2})^2
  +  (8\mu_1^2 - 2\mu_{2}) \mu_{2-r} - 4 \mu_1 \mu_{3-r} + \mu_{4-r}  \bigr ) }
}

\value{A scalar with the value of the variance estimate for the test statistic.}

\references{ Economou et. al. (2021). Hypothesis testing for the population mean and variance based on r-size biased samples, under review.}

\author{
Polychronis Economou

R implementation and documentation: Polychronis Economou  <peconom@upatras.gr>
  }


\seealso{ \code{\link{zeta_plug_in}}  }

\examples{
#s11 for the Gamma distribution for true parms=(2,3), r=1:
s11.s22(c(2,3),1, "s11", "gamma")
#s22 for for the Weibull distribution for true parms=(2,3), r=1:
s11.s22(c(2,3),1, "s22",  "weib")
}
