

\name{emxVARModel}
\alias{emxVARModel}
\alias{emxModelVAR}

\title{Create a vector autoregressive (VAR) model}

\description{
This function creates a vector autoregressive (VAR) model as an MxModel object.
}

\usage{
emxVARModel(model, data, name, run=FALSE, use, ID)
emxModelVAR(model, data, name, run=FALSE, use, ID)
}


\arguments{
   \item{model}{Currently ignored, but later will specify particular kinds of VAR models}
   \item{data}{data used for the model}
   \item{name}{character.  Optional name of the model created.}
   \item{run}{logical. Whether to run the model before returning.}
   \item{use}{character vector.  The names of the variables to use.  Currently ignored.}
   \item{ID}{character.  Name of variable that identifies each unique person.}
}

\details{
The purpose of this function is to quickly specify a vector autoregressive model.  It is currently in the early stages of development and might change considerable with regard to the \code{model} argument and the \code{ID} argument.

% The \code{data} argument can be a list of \code{data.frame} objects with one element for each \code{ID}, or a single \code{data.frame} with an \code{ID} variable that separates groups.

% When there are multiple IDs
}

\value{
    An MxModel.
}

\seealso{
\link{emxStateSpaceMixtureModel} , \link{emxFactorModel}
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
}
