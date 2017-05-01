

\name{emxResiduals}
\alias{emxResiduals}

\title{Create a residual variances matrix}

\description{
This function creates a factor loadings matrix as an MxMatrix or MxPath object.
}

\usage{
emxResiduals(x, values=.2, free=TRUE, lbound=NA, ubound=NA, path=FALSE, type='unique')
}


\arguments{
   \item{x}{character vector. The names of the variables for which residual variances are created.}
   \item{values}{numeric vector. The starting values for the variances.}
   \item{free}{logical vector. Whether the variances are free.}
   \item{lbound}{numeric vector. Lower bounds for the variances.}
   \item{ubound}{numeric vector. Upper bounds for the variances.}
   \item{path}{logical. Whether to return the MxPath object instead of the MxMatrix.}
   \item{type}{character.  The kind of residual variance structure to create.  See Details.}
}

\details{
Possible values for the \code{type} argument are 'unique' and 'identical'.  When \code{type='unique'}, each residual variances is a unique free parameter.  When \code{type='identical'}, all of the residual variances are given by a single free parameter.  In this case, all the residual variances are constrained to be equal.  However, no linear or non-liniear contraint function is used.  Rather, a single parameter occurs in multiple locations by using the same label.
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
   
# Create a residual variance matrix
require(EasyMx)
manVars <- paste0('x', 1:6)
emxResiduals(manVars, lbound=1e-6)
emxResiduals(manVars, type='identical')
emxResiduals(manVars, path=TRUE)

}
