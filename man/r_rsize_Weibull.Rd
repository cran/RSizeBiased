\name{r_rsize_Weibull}
\alias{r_rsize_Weibull}


\title{Weibull size biased random number generation of order \eqn{r} (modified).}

\description{Provides a random sample of size \eqn{n} from the \eqn{r-}size biased Weibull distribution (modified).}

\usage{r_rsize_Weibull(n,TRpar,r) }

\arguments{
  \item{n }{Number of th sample data points to be provided.}
  \item{TRpar }{A vector of length 2, containing the shape and scale parameters of the distribution.}
  \item{r }{The size (order) of the distribution. The special cases \eqn{r=1,2,3} correspond to  length, area, volume biased samples respectively and are the most frequently encountered in practice. The case \eqn{r=0} corresponds to random samples from the Weibull distribution.}
}

\details{The \eqn{r-}size random number generator from the Weibull distribution is implemented based on a change-of-variable technique, to the standard gamma distribution as described by Gove and Patil (1998).
}

\value{A vector of length \eqn{n} with the random sample.}

\references{
Gove J.H. and Patil G.P. (1998). Modeling the Basal Area-size Distribution of Forest Stands: A Compatible
Approach. Forest Science, 44(2), 285-297.
    }

\author{
Polychronis Economou

R implementation and documentation: Polychronis Economou  <peconom@upatras.gr>
  }


\seealso{ \code{\link{d_rsize_Weibull}}, \code{\link{p_rsize_Weibull}}    }

\examples{
#Random number geenration for the r-size Weibull distribution.
r_rsize_Weibull(100,c(2,3),1)
}
