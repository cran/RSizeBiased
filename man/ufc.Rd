\name{ufc}
\alias{ufc}


\title{Upper Flat Creek forest cruise tree data}

\description{Forest measurement data from the Upper Flat Creek unit of the University of Idaho Experimental Forest, measured in 1991.}

\usage{  ufc }

\format{
A data frame with 336 observations on the following 5 variables;  plot (plot label), tree (tree label),   species (species kbd with levels DF, GF, WC, WL),   dbh.cm (tree diameter at 1.37 m. from the ground, measured in centimetres.), height.m  (tree height measured in metres).
}

\details{
The inventory was based on variable radius plots with 6.43 sq. m. per ha. BAF (Basal Area Factor). The forest stand was 121.5 ha. This version of the data omits errors, trees with missing heights, and uncommon species. The four species are Douglas-fir, grand fir, western red cedar, and western larch.
}

\source{
Harold Osborne and Ross Appelgren of the University of Idaho Experimental Forest.
}





\references{
  \href{https://www.springer.com/gp/book/9781441977618}{Robinson, A.P., and J.D. Hamann. 2010. Forest Analytics with R: an Introduction. Springer.}
    }




\examples{
  data(ufc)
  }


