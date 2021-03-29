\name{d_rsize_Weibull}
\alias{d_rsize_Weibull}


\title{Weibull size biased distribution of order \eqn{r}.}

\description{Calculates the density of the \eqn{r-}size biased Weibull distribution.}

\usage{d_rsize_Weibull(x,TRpar,r) }

\arguments{
  \item{x }{Grid points where the functional is being calculated.}
  \item{TRpar }{A vector of length 2, containing the shape and scale parameters of the distribution.}
  \item{r }{The size (order) of the distribution. The special cases \eqn{r=1,2,3} correspond to  length, area, volume biased samples respectively and are the most frequently encountered in practice. The case \eqn{r=0} corresponds to random samples from the Weibull distribution.}
}

\details{The \eqn{r-}size density of the observed biased sample \eqn{X_1, \dots, X_n} is defined  by

\deqn{f_r(x; \theta)=\frac{x^r f(x; \theta)}{E(X^r)}}

where \eqn{f(x; \theta)} is the density of the Weibull distribution and \eqn{\theta} the vector of the shape and scale parameters of the distribution.
}

\value{A vector of length equal to the length of \eqn{x}.}

%\references{
%
%    }

\author{
Polychronis Economou

R implementation and documentation: Polychronis Economou  <peconom@upatras.gr>
  }


\references{ Economou et. al. (2021). Hypothesis testing for the population mean and variance based on r-size biased samples, under review.}

\seealso{ \code{\link{p_rsize_Weibull}}, \code{\link{r_rsize_Weibull}}    }

\examples{
# example of r-size Weibull distribution, r=0,1,2
x<- seq(0, 10, length=50)
dens.0.size<-d_rsize_Weibull(x,c(2,3),0)
dens.1.size<-d_rsize_Weibull(x,c(2,3),1)
dens.2.size<-d_rsize_Weibull(x,c(2,3),2)
plot(x, dens.0.size, type="l", ylab="r-denisty")
lines(x, dens.1.size, col=2)
lines(x, dens.2.size, col=3)
legend("topright",legend=c("r= 0","r= 1","r= 2"),
       col=c("black","red","green"),lty=c(1,1,1))
}
