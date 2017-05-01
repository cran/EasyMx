

\name{emxMeans}
\alias{emxMeans}

\title{Create a set of means}

\description{
This function creates a means matrix as an MxMatrix or MxPath object.
}

\usage{
emxMeans(x, values=0, free=TRUE, path=FALSE, type='saturated', name, column=TRUE, labels)
}


\arguments{
   \item{x}{character vector. The names of the variables for which covariances are created.}
   \item{values}{numeric vector. See Details.}
   \item{free}{logical vector. See Details.}
   \item{path}{logical. Whether to return the MxPath object instead of the MxMatrix.}
   \item{type}{character.  The kind of covariance structure to create.  See Details.}
   \item{name}{The name of the matrix created.}
   \item{column}{logical. Whether to create the means vector as a column or row.}
   \item{labels}{character vector. Optional labels for the means.}
}

\details{
Possible values for the \code{type} argument are 'saturated', 'equal', 'twin', 'special'.
% TODO describe the behavior for each type.
}

\value{
    Depending on the value of the \code{path} argument, either an MxMatrix or and MxPath object that can be inspected, modified, and/or included in MxModel objects.
}

\seealso{
\link{emxFactorModel}, \link{emxGrowthModel}
}

%\references{
%
%}

\examples{
   
# Create a covariance matrix
require(EasyMx)
manVars <- paste0('x', 1:6)
emxMeans(manVars, type='saturated')


}
