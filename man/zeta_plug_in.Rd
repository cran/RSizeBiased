\name{zeta_plug_in}
\alias{zeta_plug_in}


\title{\eqn{\zeta_{n,r}^i, i=1,2} test statistic for the Weibull or the gamma distribution (depending on user input.}

\description{Studentized version of the \eqn{T^i_{n,r}, i=1,2} test statistic for the Weibull/gamma distribution.}



\usage{zeta_plug_in(null_value, datain,r,EST_par,type, dist)}

\arguments{
  \item{null_value }{The parameter value in the hypothesis test under the null}
  \item{datain }{The available sample points.}
  \item{r }{The size (order) of the distribution. The special cases \eqn{r=1,2,3} correspond to  length, area, volume biased samples respectively and are the most frequently encountered in practice. The case \eqn{r=0} corresponds to random samples from the underlying distribution.}
  \item{EST_par }{A vector of length 2, containing the shape and scale parameters of the Weibull distribution.}
  \item{type }{Numeric switch: type =1 returns the \eqn{\zeta_{n,r}^1} test statistic, any other value returns \eqn{\zeta_{n,r}^2}}
  \item{dist }{ Character switch, enables the choice of distribution: type "weib" for the Weibull or "gamma" for the gamma distribution.}
}



\details{When \code{type}=1  the function returns

\deqn{\sqrt{n} \frac{T_{n,r^1} - \mu^0}{ \sigma_{1,r}(\hat \theta_n)} \rightarrow N(0,1)}

after using the fact that under the null we have \eqn{\mu_1=\mu^0}. Any other value for \code{type} returns

\deqn{\sqrt{n} \frac{T_{n,r^2} - \sigma_0^2}{ \sigma_{2,r}(\hat \theta_n)} \rightarrow N(0,1)}

in which case the fact that var(X)\eqn{=\sigma_0^2} under the null has been used.
}

\value{A scalar with the value of the test statistic.}

\references{ Economou et. al. (2021). Hypothesis testing for the population mean and variance based on r-size biased samples, under review.}

\author{
Polychronis Economou

R implementation and documentation: Polychronis Economou  <peconom@upatras.gr>
  }


%\seealso{ \code{\link{zeta_plug_ing}}  }

\examples{
data(ufc)
datain_r <- ufc[,4]
nullMEAN <- 14
# ml estimates = c(2.6555,8.0376),  taken from section 6.2 in Economou et. al. (2021).
zeta_plug_in(nullMEAN, datain_r, 2, c(2.6555,8.0376),1, "gamma") #corresponds to mean

nullVar <- 180
zeta_plug_in(nullVar, datain_r, 2, c(2.6555,8.0376),2, "gamma") #corresponds to var

}
