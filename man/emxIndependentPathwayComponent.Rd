\name{emxIndependentPathwayComponent}
\alias{emxIndependentPathwayComponent}

\title{Creates component for a Biometric Independent Pathway Model}

\description{
UNDER ACTIVE DEVELOPMENT. DO NOT TRUST. This function creates all the objects needed for an Independent Pathway component.
}

\usage{
emxIndependentPathwayComponent(x, xname, xvalues=.8, xfree=TRUE, xlbound=NA, xubound=NA,
                               h=2, hname=paste0('H', xname), hvalues, hlabels)
}

\arguments{
   \item{x}{character vector. The base names of the variables used for the matrix with no repetition for twins (x, y, z not x1, y1, z1, x2, y2, z2).}
   \item{xname}{character. Name of the component matrix.}
   \item{xvalues}{Matrix of fixed and/or starting values of parameters}
   \item{xfree}{Matrix of TRUE (for free) and FALSE (for fixed) values}
   \item{xlbound}{Matrix of numeric lower bounds for parameters}
   \item{xubound}{Matrix of numeric upper bounds for parameters}
   \item{h}{numeric. The number of variables for the relatedness matrix, i.e. the number of critters with relationships}
   \item{hname}{character. Name of the relatedness matrix.}
   \item{hvalues}{numeric vector. Values for the relatedness matrix.}
   \item{hlabels}{character vector. Labels for the relatedness matrix.}
}

\details{
This function is a combination of emxCholeskyVariance, emxRelatednessMatrix, and emxKroneckerVariance.
}

\value{
    A list with elements (1) the lower triangular matrix for the Cholesky, (2) the full positive definite variance matrix, (3) the relatedness matrix, and (4) the Kronecker product of the variance matrix and the relatedness matrix.
}

\seealso{
\link{emxGeneticFactorComponent}
}

%\references{
%
%}

\examples{
   
# Create an ACE model in 22 lines
require(EasyMx)
require(OpenMx)
data(twinData)
twinVar = names(twinData)
selVars <- c('ht1', 'bmi1','ht2','bmi2')
mzdzData <- subset(twinData, zyg \%in\% c(1, 3), c(selVars, 'zyg'))
mzdzData$RCoef <- c(1, NA, .5)[mzdzData$zyg]
nVar = length(selVars)/2
x <- paste0('x', 1:nVar)

acomp <- emxCholeskyComponent(x, 'A', hvalues=c(1, .5, 1), hlabels=c(NA, 'data.RCoef', NA))
ccomp <- emxCholeskyComponent(x, 'C', hvalues=c(1, 1, 1))
ecomp <- emxCholeskyComponent(x, 'E', hvalues=c(1, 0, 1))
totalVar <- mxAlgebra(AKron + CKron + EKron, 'V', dimnames=list(selVars, selVars))
totalMean <- emxMeans(selVars, type='twin')
expect <- mxExpectationNormal(totalVar$name, totalMean$name)
fitfun <- mxFitFunctionML()

comlist <- c(acomp, ccomp, ecomp, list(totalVar, totalMean, expect, fitfun))

model <- mxModel('model', comlist, mxData(mzdzData, 'raw'))
\dontrun{
run2 <- mxRun(model)
}

}
