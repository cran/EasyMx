\name{emxGeneticFactorComponent}
\alias{emxGeneticFactorComponent}

\title{Creates component for a Genetic Factor Model}

\description{
This function creates all the objects needed for Genetic Factor component
}

\usage{
emxGeneticFactorComponent(x, xname, xvalues=.8, xfree=TRUE, xlbound=NA, xubound=NA,
                          h=2, hname=paste0('H', xname), hvalues, hlabels)
}

\arguments{
   \item{x}{character vector. The base names of the variables used for the matrix with no repetition for twins (x, y, z not x1, y1, z1, x2, y2, z2).}
   \item{xname}{character. Name of the component matrix.}
   \item{xvalues}{numeric vector. Values of the genetic factor loadings.}
   \item{xfree}{logical vector. Whether the genetic factor loadings are free.}
   \item{xlbound}{numeric vector. Lower bounds of the factor loadings.}
   \item{xubound}{numeric vector. Upper bounds of the factor loadings.}
   \item{h}{numeric. The number of variables for the relatedness matrix, i.e. the number of critters with relationships}
   \item{hname}{character. Name of the relatedness matrix.}
   \item{hvalues}{numeric vector. Values for the relatedness matrix.}
   \item{hlabels}{character vector. Labels for the relatedness matrix.}
}
\details{
This function is a combination of emxGeneticFactorVariance, emxRelatednessMatrix, and emxKroneckerVariance.
}

\value{
    A list with elements (1) the genetic factor loadings matrix, (2) the full positive definite variance matrix, (3) the relatedness matrix, and (4) the Kronecker product of the variance matrix and the relatedness matrix.
}

\seealso{
\link{emxCholeskyComponent}
}

%\references{
%
%}

\examples{
   
# Create genetic factor A component for DZ twins
require(EasyMx)
xvars <- paste0('x', 1:4)
acomp <- emxGeneticFactorComponent(xvars, 'A', hvalues=c(1, .5, 1))

}
