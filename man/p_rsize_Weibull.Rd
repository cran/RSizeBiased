\name{p_rsize_Weibull}
\alias{p_rsize_Weibull}


\title{Weibull size biased c.d.f. of order \eqn{r}.}

\description{Calculates the cumulative distribution of the \eqn{r-}size biased Weibull distribution.}

\usage{p_rsize_Weibull(q,TRpar,r) }

\arguments{
  \item{q }{Points where the functional is being calculated.}
  \item{TRpar }{A vector of length 2, containing the shape and scale parameters of the distribution.}
  \item{r }{The size (order) of the distribution. The special cases \eqn{r=1,2,3} correspond to  length, area, volume biased samples respectively and are the most frequently encountered in practice. The case \eqn{r=0} corresponds to random samples from the Weibull distribution.}
}

\details{The \eqn{r-}size c.d.f. of the Weibull density is defined  by

\deqn{F_r(y; \theta)=\int_{0}^{y} \frac{x^r f(x; \theta)}{E(X^r)} \,dx}

where \eqn{\theta} is a bivariate vector with the the shape and scale of the Weibull distribution.
}

\value{A vector of length equal to the lemgth of \eqn{x}.}

%\references{
%
%    }

\author{
Polychronis Economou

R implementation and documentation: Polychronis Economou  <peconom@upatras.gr>
  }


\references{ Economou et. al. (2021). Hypothesis testing for the population mean and variance based on r-size biased samples, under review.}

\seealso{ \code{\link{d_rsize_Weibull}}, \code{\link{r_rsize_Weibull}}    }

\examples{
# c.d.f of the r-size Weibull distribution, r=0,1,2 evalutated at a specific point x.
x<- 2
dist.0.size<-p_rsize_Weibull(x,c(2,3),0)
dist.1.size<-p_rsize_Weibull(x,c(2,3),1)
dist.2.size<-p_rsize_Weibull(x,c(2,3),2)
}
