

\name{emxRelatednessMatrix}
\alias{emxRelatednessMatrix}

\title{Create a relatedness matrix}

\description{
This function creates a relatedness matrix as an MxMatrix, often used in behavior genetics modeling.
}

\usage{
emxRelatednessMatrix(nvar, values, labels, name='h')
}


\arguments{
   \item{nvar}{numeric. The number of variables for the matrix, i.e. the number of rows or columns.}
   \item{values}{numeric vector. Values used in the matrix.}
   \item{labels}{character vector. Labels of the elements in the matrix.  See Details.}
   \item{name}{character. The name of the matrix created.}
}

\details{
The \code{labels} argument can be used to create a "definition variable" which populates the value from one of the data columns for each row in the data.  In this context, if the genetic relatedness coefficient between a pair of individuals is given by a column in the data then that information can be used to create in the relatedness matrix.  Alternatively, multiple groups can be created 
}

\value{
    An MxMatrix object.
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
ahmat <- emxRelatednessMatrix(2, c(1, .5, 1), labels=c(NA, 'data.RCoef', NA), name='AH')
# data.RCoef creates a definition variable and ignores the .5 value.

}
