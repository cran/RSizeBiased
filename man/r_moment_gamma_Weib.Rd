\name{r_moment_gamma_Weib}
\alias{r_moment_gamma_Weib}


\title{\eqn{r-}th moment of the gamma or the Weibull distribution.}

\description{Calculates the \eqn{r-}th moment of the gamma or Weibull distribution.}

\usage{r_moment_gamma_Weib(TRpar,r,dist)}

\arguments{
\item{TRpar }{A vector of length 2, containing the shape and scale parameters of the distribution.}
\item{r }{The size (order) of the distribution. The special cases \eqn{r=1,2,3} correspond to  length, area, volume biased samples respectively and are the most frequently encountered in practice. The case \eqn{r=0} corresponds to random samples from the Gamma distribution.}
 \item{dist }{ Character switch, enables the choice of distribution: type "weib" for the Weibull or "gamma" for the gamma distribution.}
}

\details{In the case of the \eqn{\Gamma(\alpha, \beta)} distribution the \eqn{r}-th moment is given by

\deqn{\mu_r = \int_0^{\infty} x^r f(x;\alpha, \beta)\,dx =\beta^r \frac{\Gamma(\alpha+r)}{\Gamma(\alpha)},  \alpha> -r }

while for the \eqn{W(\alpha, \beta)} distribution the \eqn{r}-th moment is given by

\deqn{\mu_r = \int_0^{\infty} x^r f(x;\alpha, \beta)\,dx = \beta^r \Gamma\left(1+\frac{\alpha}{r}\right), \alpha> -r }



}

\value{A scalar, the value of the moment.}

\references{ Economou et. al. (2021). Hypothesis testing for the population mean and variance based on r-size biased samples, under review.}

\author{
Polychronis Economou

R implementation and documentation: Polychronis Economou  <peconom@upatras.gr>
  }




\examples{
#r-moment for the Gamma distribution for true parms=(2,3), r=1:
r_moment_gamma_Weib(c(2,3),1, "gamma")
#r-moment for for the Weibull distribution for true parms=(2,3), r=1:
r_moment_gamma_Weib(c(2,3),1, "weib")
}
