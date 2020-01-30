

\name{emxFactorModel}
\alias{emxFactorModel}
\alias{emxModelFactor}

\title{Create a factor model}

\description{
This function creates a factor model as an MxModel object.
}

\usage{
emxFactorModel(model, data, name, run=FALSE, identification, use, ordinal,
    ..., parameterization=c("lisrel", "ifa"), weight = as.character(NA))
emxModelFactor(model, data, name, run=FALSE, identification, use, ordinal,
    ..., parameterization=c("lisrel", "ifa"), weight = as.character(NA))
}


\arguments{
   \item{model}{named list. Gives the factor loading pattern. See Details.}
   \item{data}{data used for the model}
   \item{name}{character.  Optional name of the model created.}
   \item{run}{logical. Whether to run the model before returning.}
   \item{identification}{Not yet implemented.  How the model is identified.  Currently ignored.}
   \item{use}{character vector.  The names of the variables to use.}
   \item{ordinal}{character vector.  The names of the ordinal variables.}
   \item{...}{Force later arguments to be named.  ... is ignored.}
   \item{parameterization}{character.  Whether to specify the model as a LISREL SEM or as an Item Factor Analysis}
   \item{weight}{character. Name of the data column used for sample weights.}
}

\details{
The \code{model} argument must be a named list.  The names of the list give the names of the latent variables.  Each list element gives the names of the variables that load onto that latent variable.  This may sound complicated, but the example below makes this more clear.  It is intended to be visually intuitive.


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
xmap <- list(F1=paste0('x', 1:6), F2=paste0('y', 1:3), F3=paste0('z', 1:3))
\dontrun{
mod <- emxFactorModel(xmap, data=myFADataRaw, run=TRUE)
}

}
