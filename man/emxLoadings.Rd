

\name{emxLoadings}
\alias{emxLoadings}

\title{Create a factor loadings matrix}

\description{
This function creates a factor loadings matrix as an MxMatrix or MxPath object.
}

\usage{
emxLoadings(x, values=.8, free=TRUE, path=FALSE)
}


\arguments{
   \item{x}{named list. Gives the factor loading pattern. See Details.}
   \item{values}{numeric vector. The starting values for the nonzero loadings.}
   \item{free}{logical vector. Whether the nonzero loadings are free.}
   \item{path}{logical. Whether to return the MxPath object instead of the MxMatrix.}
}

\details{
The \code{x} argument must be a named list.  The names of the list give the names of the latent variables.  Each list element gives the names of the variables that load onto that latent variable.  This may sound complicated, but the example below makes this more clear.  It is intended to be visually intuitive.
}

\value{
    Depending on the value of the \code{path} argument, either an MxMatrix or and MxPath object that can be inspected, modified, and/or included in MxModel objects.
}

\seealso{
\link{emxFactorModel}
}

%\references{
%
%}

\examples{
   
# Create a loadings matrix
require(EasyMx)
xmap <- list(F1=paste0('x', 1:6), F2=paste0('y', 1:3), F3=paste0('z', 1:3))
emxLoadings(xmap)
emxLoadings(xmap, path=TRUE)

}
