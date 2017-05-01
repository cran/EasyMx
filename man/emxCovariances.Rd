

\name{emxCovariances}
\alias{emxCovariances}

\title{Create a set of covariances}

\description{
This function creates a covariance matrix as an MxMatrix or MxPath object.
}

\usage{
emxCovariances(x, values, free, path=FALSE, type, name='Variances')
}


\arguments{
   \item{x}{character vector. The names of the variables for which covariances are created.}
   \item{values}{numeric vector. See Details.}
   \item{free}{logical vector. See Details.}
   \item{path}{logical. Whether to return the MxPath object instead of the MxMatrix.}
   \item{type}{character.  The kind of covariance structure to create.  See Details.}
   \item{name}{The name of the matrix created.}
}

\details{
Possible values for the \code{type} argument are 'independent', 'full', and 'corr'.  When \code{type='independent'}, the remaining arguments are passes to \code{\link{emxResiduals}}.  The \code{values} and \code{free} arguments are only used when the \code{type} argument is 'independent'.  For all other cases, they are ignored.

When \code{type='full'}, a full covariance matrix is created.  That is, a symmetric matrix is created with all unique elements freely estimated.  The starting values for the variances are all 1; for the covariances, all 0.5.

When \code{type='corr'}, a full correlation matrix is created.  That is, a symmetric matrix is created with all unique elements not on the diagonal freely estimated.  The starting values for the correlations are all 0.5.  The variances are fixed at 1.

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
latVars <- paste0('F', 1:2)
emxCovariances(manVars, type='full')
emxCovariances(latVars, type='corr', path=TRUE)


}
