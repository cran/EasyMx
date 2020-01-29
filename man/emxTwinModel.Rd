\name{emxTwinModel}
\alias{emxTwinModel}
\alias{emxModelTwin}

\title{Creates behavior genetics Twin Model}

\description{
This function creates an MxModel and associated objects for a basic Twin model.
}

\usage{
emxTwinModel(model, relatedness, data, run=FALSE, use, name='model', components='ACE')
emxModelTwin(model, relatedness, data, run=FALSE, use, name='model', components='ACE')
}


\arguments{
   \item{model}{Description of the model. Currently ignored.}
   \item{relatedness}{Description of the relatedness patterns.  Currently the name of the variable that gives the coefficient of relatedness.}
   \item{data}{data.frame or matrix. The data set used in the model.}
   \item{run}{logical. Whether to run the model before returning.}
   \item{use}{character vector. Names of the variables used in the model.}
   \item{name}{character. Name of the model.}
   \item{components}{character. Name of the variance components to include. Current valid options are 'ACE' and 'ADE'}
}

\details{
Because the model argument is ignored and the relatedness argument has limited use, this function only constructs a very basic and rigid Twin model.  It creates a Cholesky model with A, C, and E components or a Cholesky model with A, D, and E components. The means are constrained equal across twins.
}

\value{
    MxModel.
}

\seealso{
\link{emxFactorModel}
}

%\references{
%
%}

\examples{
   
# Create an ACE model in 10 lines
# 8 of those are data handling.
# 2 are the actual model.
require(EasyMx)
require(OpenMx)
data(twinData)
twinVar = names(twinData)
selVars <- c('ht1', 'bmi1','ht2','bmi2')
mzdzData <- subset(twinData, zyg \%in\% c(1, 3), c(selVars, 'zyg'))
mzdzData$RCoef <- c(1, NA, .5)[mzdzData$zyg]

run3 <- emxTwinModel(model='Cholesky', relatedness='RCoef',
	data=mzdzData, use=selVars, run=TRUE, name='TwCh')


}
