\name{emxKroneckerVariance}
\alias{emxKroneckerVariance}

\title{Creates a large Variance matrix by Kroneckering two smaller matrices}

\description{
This function creates the wide format variance matrix when combined with a relatedness matrix
}

\usage{
emxKroneckerVariance(h, v, name)
}


\arguments{
   \item{h}{MxMatrix. Left hand side of the Kronecker product.  Typically the relatedness matrix.}
   \item{v}{MxMatrix. Right hand side of the Kronecker product.  Typically the variance matrix.}
   \item{name}{character. Name of the resulting large matrix.}
}

\details{
In many behavior genetic models, a relationship matrix is combined with a base variance matrix.  The combination is done with a Kronecker product so that the variance exists (possibly weighted by zero or another number) for each member of the relationship.
}

\value{
    MxAlgebra
}

\seealso{
\link{emxRelatednessMatrix}
}

%\references{
%
%}

\examples{
   
# Create a loadings matrix
require(EasyMx)
x <- letters[23:26]
amat <- emxCholeskyVariance(x, 'A')
ahmat <- emxRelatednessMatrix(2, c(1, .5, 1), name='AH')
ab <- emxKroneckerVariance(ahmat, amat[[2]], 'AB')

}
