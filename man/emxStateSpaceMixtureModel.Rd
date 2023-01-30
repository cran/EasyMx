

\name{emxStateSpaceMixtureModel}
\alias{emxStateSpaceMixtureModel}
\alias{emxModelStateSpaceMixture}

\title{Create a state space mixture model}

\description{
This function creates a state space mixture model as an MxModel object.
}

\usage{
emxStateSpaceMixtureModel(model, data, name, run=FALSE, use, ID, time, ...)
emxModelStateSpaceMixture(model, data, name, run=FALSE, use, ID, time, ...)
}


\arguments{
   \item{model}{list of MxModel objects, each of which should be a state space model}
   \item{data}{data used for the model}
   \item{name}{character.  Optional name of the model created.}
   \item{run}{logical. Whether to run the model before returning.}
   \item{use}{character vector.  The names of the variables to use.  Currently ignored.}
   \item{ID}{character.  Name of variable that identifies each unique person.}
   \item{time}{character.  Name of the variable that gives the time of each obsevation.  Currently ignored.}
   \item{...}{Force later arguments to be named.  ... is ignored.}
}

\details{
The idea of state space mixture modeling is to model multiple, multivariate time series while allowing for qualitative differences between the time series.  Suppose you have a multivariate time series for several people.  You think some people should have the same time series model, but not everyone.  You think there should be a small number of homogeneous sets of people that follow the same time series model, but you do not know which people or the exact parameter values of the candidate time series models.  This function presents one solution to this problem.

State space mixture modeling begins with a set of candidate state space models, and uses these state space models as the mixture classes.  The goal is to simultaneously estimate the free parameters of the state space models, and estimate which multivariate time series (e.g., perseon) belongs to which mixture class.

The component state space models may share some free parameters or none.  Note that free parameters with the same name are constrained to be equal across all models.  Conversely, unnamed free parameters are given unique names and are allowed to differ for each person-mixture combination, which creates a very large number of free parameters.  We strongly encourage you to name all of your free parameters in the \code{model} list to avoid melting your computer's CPU.

The \code{model} argument currently must be a list.  The elements of the list should be MxModel objects.  Each list element forms a mixture class in the final model.

This function creates a multigroup mixture model where the mixture classes are the elements of the \code{model} list argument.  Each unique \code{ID} forms a group.

The \code{data} argument can be a list of \code{data.frame} objects with one element for each \code{ID}, or a single \code{data.frame} with an \code{ID} variable that separates groups.
}

\value{
    An MxModel.
}

\seealso{
\link{emxStateSpaceMixtureClassify} , \link{emxMixtureModel}
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
# Make and fit the state space mixture model
ssmm <- emxStateSpaceMixtureModel(model=list(vm1, vm2),
    data=ds1, ID='id', run=TRUE)
}
}
