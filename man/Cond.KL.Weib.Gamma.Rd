\name{Cond.KL.Weib.Gamma}
\alias{Cond.KL.Weib.Gamma}


\title{Kullback-Leibler divergence between the (parametrized with respect to shape and mean or variance) of the Weibull or gamma distribution and its (assumed) maximum likelihood estimates.}

\description{The function returns the Kullback-Leibler divergence (minus a constant) between the (parametrized with respect to shape and mean or variance) underlying Weibull or gamma distribution and its (assumed) maximum likelihood estimates.}

\usage{Cond.KL.Weib.Gamma(par,nullvalue,hata,hatb,type,dist)}


\arguments{
  \item{par }{The (actual) shape parameter \eqn{\alpha} of the distribution.}
  \item{nullvalue }{The (actual) distribution mean or variance.}
   \item{hata }{Maximum likelihood estimate of the shape parameter of the distribution.}
  \item{hatb }{Maximum likelihood estimate of the scale parameter of the distribution.}
   \item{type }{Numeric switch, enables the choice of mean or variance: type: 1 for mean, 2 (or any other value != 1) for variance.}
   \item{dist }{ Character switch, enables the choice of distribution: type "weib" for the Weibull or "gamma" for the gamma distribution.}
}

\details{The Kullback-Leibler divergence between the Weibull(\eqn{\alpha, \beta)} or the gamma(\eqn{\alpha, \beta)} and its maximum likelihood estimate Gamma(\eqn{\hat \alpha, \hat \beta)} is  given by

\deqn{
D_{KL} =  (\hat \alpha -1)\Psi(\hat \alpha) - \log\hat \beta - \hat \alpha - \log \Gamma(\hat \alpha) + \log\Gamma(  \alpha) + \alpha \log \beta - (\alpha -1)(\Psi(\hat \alpha)  + \log \hat \beta) + \frac{ \hat \beta \hat \alpha}{\lambda}.
}
Since \eqn{D_{KL}} is  used to determine the closest distribution - given its mean or variance - to the estimated gamma p.d.f., the first four terms are omitted from the function outcome, i.e. the function returns the result of the following quantity:

\deqn{
 \log\Gamma(  \alpha) + \alpha \log \beta - (\alpha -1)(\Psi(\hat \alpha)  + \log \hat \beta) + \frac{ \hat \beta \hat \alpha}{\lambda}.
}

For the Weibull distribution the corresponding formulas are

\deqn{
D_{KL} =   \log \frac{\hat \alpha}{{\hat \beta}^{\hat \alpha}} - \log \frac{\alpha}{{\beta}^{\alpha}} + (\hat \alpha - \alpha) \left ( \log \hat \beta - \frac{\gamma}{\hat \alpha} \right ) + \left (\frac{\hat \beta}{\beta} \right )^\alpha \Gamma\left ( \frac{\alpha}{\hat \alpha} +1 \right ) -1
}
and since \eqn{D_{KL}} is  used to determine the closest distribution - given its mean or variance - to the estimated gamma p.d.f., the first   term is omitted from the function outcome, i.e. the function returns the result of the following quantity:

\deqn{
   - \log \frac{\alpha}{{\beta}^{\alpha}} + (\hat \alpha - \alpha) \left ( \log \hat \beta - \frac{\gamma}{\hat \alpha} \right ) + \left (\frac{\hat \beta}{\beta} \right )^\alpha \Gamma\left ( \frac{\alpha}{\hat \alpha} +1 \right ) -1
}


}

\value{A scalar, the value of the Kullback-Leibler divergence (minus a constant).}

\references{ Economou et. al. (2021). Hypothesis testing for the population mean and variance based on r-size biased samples, under review.}

\author{
Polychronis Economou

R implementation and documentation: Polychronis Economou  <peconom@upatras.gr>
  }


%\seealso{ \code{\link{KLGamma_for_min_given_variance}}     }

\examples{
#K-L divergence for the Gamma distribution for shape=2
#and variance=3 and their assumed MLE=(1,1):
 Cond.KL.Weib.Gamma(2,3,1,1,2, "gamma")
#K-L divergence for the Weibull distribution for shape=2
#and variance=3 and their assumed MLE=(1,1):
 Cond.KL.Weib.Gamma(2,3,1,1,2, "weib")
}%
%}
