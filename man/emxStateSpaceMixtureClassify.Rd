

\name{emxStateSpaceMixtureClassify}
\alias{emxStateSpaceMixtureClassify}

\title{Classify time series in a state space mixture model}

\description{
This function classifies time series (usually people) in a state space mixture model.
}

\usage{
emxStateSpaceMixtureClassify(model)
}


\arguments{
   \item{model}{MxModel.  The output from \code{\link{emxStateSpaceMixtureModel}}}
}

\details{
This is a helper function for state space mixture modeling.  The function will almost exclusively be used in conjunction with emxStateSpaceMixtureModel.  The present function takes a state space mixture model as input, and returns detailed information about the most likely class for each unique ID.
}

\value{
    A named list with elements
    \item{estimated_classes}{A vector of the most likely class for each person.  Dimension is people.}
    \item{joint_m2ll}{A matrix of joint minus two summed log likelihoods of each person \emph{and} each class.  Dimension is people by classes.}
    \item{m2ll}{A matrix of minus two summed log likelihoods of each person \emph{given} each class.  Dimension is people by classes.}
    \item{likelihood}{An array of the likelihoods (i.e., probability densities) of each combination of time point, person, and class.  Dimension is time points by people by classes.}
}

\seealso{
\link{emxStateSpaceMixtureModel} , \link{emxMixtureModel}
}
% \link{emxStateSpaceModel}

%\references{
%
%}

\examples{
   
# Example
require(EasyMx)
data(myFADataRaw)

ds0 <- myFADataRaw[,1:3]

# Make a VAR Model
vm <- emxVARModel(data=ds0, use=names(ds0), name='varmodel')

# Re-label parameters to have different AR parameters
# for class 1 and class 2
vm1 <- OpenMx::omxSetParameters(vm, labels=vm$Dynamics$labels,
    newlabels=paste0(vm$Dynamics$labels, '_k1'), name='klass1')
vm2 <- OpenMx::omxSetParameters(vm, labels=vm$Dynamics$labels,
    newlabels=paste0(vm$Dynamics$labels, '_k2'), name='klass2')

# Pretend you have a data set of 50 people
# each measured 10 times on 3 variables
ds1 <- myFADataRaw[, 1:3]
ds1$id <- rep(1:50, each=nrow(myFADataRaw)/50)

\dontrun{
# Make the state space mixture model
ssmm <- emxStateSpaceMixtureModel(model=list(vm1, vm2),
    data=ds1, ID='id')

# Fit model
ssmmr <- mxRun(ssmm)

# Extract estimated classes and diagnostics
eclasses <- emxStateSpaceMixtureClassify(ssmmr)
}

}
