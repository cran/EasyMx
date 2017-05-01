

\name{emxMixtureModel}
\alias{emxMixtureModel}
\alias{emxModelMixture}

\title{Create a mixture model}

\description{
This function creates a mxiture model as an MxModel object.
}

\usage{
emxMixtureModel(model, data, run=FALSE, p=NA, ...)
emxModelMixture(model, data, run=FALSE, p=NA, ...)
}


\arguments{
   \item{model}{list. The MxModel objects that compose the mixture.}
   \item{data}{data used for the model}
   \item{run}{logical. Whether to run the model before returning.}
   \item{p}{character.  Optional name of the mixing proportions matrix.}
   \item{...}{Further Mx Objects passed into the mixture model.}
}

\details{
The \code{model} argument is list of MxModel objects.  These are the classes over which the mixture model operates.

The \code{p} argument is optional.  If not specified, the function will create and properly scale the mixing proportions for you.  If specified, the Mx Object that gives the mixing proportions should be a column vector (one-column matrix).

}

\value{
    An MxModel.
}

\seealso{
\link{emxLoadings}
}

%\references{
%
%}

\examples{
   
# Example
require(EasyMx)
data(myFADataRaw)
xmap1 <- list(F1=paste0('x', 1:6), F2=paste0('y', 1:3), F3=paste0('z', 1:3))
mod1 <- emxFactorModel(xmap1, data=myFADataRaw, name='m1')

xmap2 <- list(F1=c(paste0('x', 1:6), paste0('y', 1:3), paste0('z', 1:3)))
mod2 <- emxFactorModel(xmap2, data=myFADataRaw, name='m2')

mod <- emxMixtureModel(list(mod1, mod2), data=myFADataRaw)
# To estimate parameters either
#  1. mod <- mxRun(mod)   or
#  2. include run=TRUE in the arguments above
summary(mod)
coef(mod)


}
