\name{Size.BiasedMV.Tests}
\alias{Size.BiasedMV.Tests}


\title{Test statistics.}

\description{The function returns the test statistics for testing a null hypothesis for the mean and a null hypothesis for the varaince.}

\usage{Size.BiasedMV.Tests(datain_r,r,nullMEAN,nullVAR,start_par,nboot,alpha,prior_sel,distr)}




\arguments{
 \item{datain_r}{The available sample points.}
\item{r}{The size (order) of the distribution. The special cases \eqn{r=1,2,3} correspond to  length, area, volume biased samples respectively and are the most frequently encountered in practice. The case \eqn{r=0} corresponds to random samples from the gamma or theWeibull distribution.}
\item{nullMEAN}{The null value of the distribution mean.}
\item{nullVAR}{The null value of the distribution variance.}
\item{start_par}{Vector with two values, containing the starting values for the MLE for the two parameter distribution (Weibull or gamma) .}
\item{nboot}{Defines the  number of bootstrap replications.}
\item{alpha}{Significance level.}
\item{prior_sel}{"normal" for the normal distribution or "gamma" for the gamma.}
\item{distr }{ Character switch, enables the choice of distribution: type "weib" for the Weibull or "gamma" for the gamma distribution.}
}


\details{The test statistics implemented are given by the Plug-in and the bootstrap Methods as described in section 3.1 and 3.2 of Economou et al (2021).
}

\value{An object containing the following components.
\item{par}{A vector of the MLE of the distribution parameters.}
\item{loglik}{A scalar, the maximized log-likelihood.}
\item{CovMatrix}{The Variance - Covariance matrix of the MLEs.}
\item{Zeta_i}{A vector of the values of the \eqn{\zeta_{n,r}^i, i=1,2} test statistics (if defined)}
\item{Tivalues}{A vector of the values of the \eqn{T^i_{n,r}, i=1,2} test statistics}
\item{T1_bootstrap_quan}{A vector of the bootstrap quantiles for the \eqn{T^1_{n,r}} test statistic for each one of the significance levels alpha.}
\item{T2_bootstrap_quan}{A vector of the bootstrap quantiles for the \eqn{T^2_{n,r}} test statistic for each one of the significance levels alpha.}
\item{NullValues}{A vector of the null values of the distribution mean and variance.}
\item{distribution}{Character representing the choice of distribution: "weib" for the Weibull or "gamma" for the gamma distribution.}
\item{alpha}{A vector of significance levels for the test level.}
\item{bootstrap_p_mean}{A scalar with the bootstrap p-value for testing the mean.}
\item{bootstrap_p_var}{A scalar with the bootstrap p-value for testing the variance.}
\item{decision}{A matrix of 0 and 1 of the decisions taken for each one of the significance levels alpha based on the bootstrap method. The first row corresponds to the null hypothesis for the mean and the second to the null hypothesis for the variance.}
\item{asymptotic_p_mean}{A scalar with the asymptotic p-value for testing the mean (if \eqn{\zeta_{n,r}^1} is defined).}
\item{asymptotic_p_var}{A scalar with the asymptotic p-value for testing the variance (if \eqn{\zeta_{n,r}^2} is defined).}
\item{decisionasympt}{A matrix of 0 and 1 of the decisions taken for each one of the significance levels alpha based on the plug-in method and the asymptotic distribution of the test statistics. The first row corresponds to the null hypothesis for the mean and the second to the null hypothesis for the variance.}
\item{prior_selection}{Character representing the choice of the prior distribution for the bootstrap method: "normal" for the normal distribution or "gamma" for the gamma.}
}

\references{Economou et. al. (2021). Hypothesis testing for the population mean and variance based on r-size biased samples, under review.}

\author{
Polychronis Economou

R implementation and documentation: Polychronis Economou  <peconom@upatras.gr>
  }


%\seealso{ \code{\link{KLWeibull_for_min_given_mean}}     }

\examples{
data(ufc)
datain_r <- ufc[,4]
nullMEAN <- 14 #according to null mean in Sec. 6.3,  Economou et. al. (2021).
nullVAR <- 180 #according to null variance in Sec. 6.3,  Economou et. al. (2021).
Size.BiasedMV.Tests(datain_r, 2, nullMEAN, nullVAR,  c(2,3), 100, 0.05, "normal", "gamma")
}
