

\name{emxCholeskyVariance}
\alias{emxCholeskyVariance}

\title{Create a variance matrix in Cholesky form}

\description{
This function creates a Cholesky variance matrix and associated MxMatrix and MxAlgebra objects.
}

\usage{
emxCholeskyVariance(x, name, values=.8, free=TRUE)
}


\arguments{
   \item{x}{character vector. The names of the variables used for the matrix.}
   \item{name}{character. The name of the variance matrix created.}
   \item{values}{numeric vector. The starting values for the lower triangular matrix.}
   \item{free}{logical vector. Whether the lower triangular elements are free.}
}

\details{
This is a helper function for creating a matrix that is symmetric and positive definite. Full covariance matrices are the most common case of these.  In a behavior genetics modeling context, Cholesky components can be created for Additive genetics, Common environments, and unique Environments. These are unrestrictive models of the covariances of multiple phenotypes.
}

\value{
    A list with two components. The first component is the lower triangular MxMatrix. The second component is an MxAlgebra, the result of which is the positive definite variance matrix.
}

\seealso{
\link{emxGeneticFactorVariance}
}

%\references{
%
%}

\examples{
   
# Create a Cholesky variance matrix called 'A'
require(EasyMx)
nVar <- 3
x <- paste0('x', 1:nVar)
amat <- emxCholeskyVariance(x, 'A')

}
