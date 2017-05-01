

\name{emxThresholds}
\alias{emxThresholds}

\title{Create a set of thresholds for ordinal data}

\description{
This function creates a threshold matrix as an MxMatrix object.
}

\usage{
emxThresholds(data, ordinalCols)
}


\arguments{
   \item{data}{The data frame or matrix for which thresholds should be created.}
   \item{ordinalCols}{optional character vector.  The names of the ordinal variables in the data.}
}

%\details{
% TODO add info about the deviations used to keep thresholds in order.
%}

\value{
    An MxMatrix giving the thresholds.
}

\seealso{
\link{emxFactorModel}, \link{emxGrowthModel}
}

%\references{
%
%}

\examples{
   
# Example
require(EasyMx)
data(jointdata)
jointdata[, c(2, 4, 5)] <- mxFactor(jointdata[,c(2, 4, 5)],
	levels=sapply(jointdata[,c(2, 4, 5)], function(x){sort(unique(x))}))
emxThresholds(jointdata, c(FALSE, TRUE, FALSE, TRUE, TRUE))

}
