\name{T1T2.Mean.Var}
\alias{T1T2.Mean.Var}


\title{Test statistic \eqn{T_{n,r}^1} or \eqn{T_{n,r}^2} depending on user input.}

\description{The test statistics \eqn{T_{n,r}^1} and \eqn{T_{n,r}^2}  are consistent estimators of the mean value \eqn{\mathrm{E}(X)} and variance   \eqn{\mathrm{Var}(X)} respectively given an \eqn{r-}size biased sample.}

\usage{T1T2.Mean.Var(datain,r, type) }

\arguments{
  \item{datain }{The available sample points.}
  \item{r }{The size (order) of the distribution. The special cases \eqn{r=1,2,3} correspond to  length, area, volume biased samples respectively and are the most frequently encountered in practice. The case \eqn{r=0} corresponds to random samples from the underlying distribution.}
  \item{type }{Numeric switch: type =1 corresponds to the T1 statistic while any other numeric value will cause calculation of T2.}
}

\details{ The test statistic \eqn{T_{n,r}^1} is defined by

\deqn{ T_{n,r}^{1}=\frac{\sum_{i=1}^n X_i^{1-r}}{\sum_{i=1}^n X_i^{-r}}.}

The test statistic \eqn{T_{n,r}^2} is defined by

\deqn{ T_{n,r}^{2}= \frac{\sum_{i=1}^n X_i^{2-r}}{\sum_{i=1}^nX_i^{-r}}-{\left(\frac{\sum_{i=1}^n X_i^{1-r}}{\sum_{i=1}^n X_i^{-r}}\right)^2}.}

}

\value{A scalar, the value of the test statistic for the given sample.}

\references{ Economou et. al. (2021). Hypothesis testing for the population mean and variance based on r-size biased samples, under review.}

\author{
Polychronis Economou

R implementation and documentation: Polychronis Economou  <peconom@upatras.gr>
  }


\examples{
#e.g.:
T1T2.Mean.Var(rgamma(100, 2,3),0, 1)
}
