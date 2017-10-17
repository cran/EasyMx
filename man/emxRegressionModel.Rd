

\name{emxRegressionModel}
\alias{emxRegressionModel}
\alias{emxModelRegression}

\title{Create a regression model}

\description{
This function creates a regression model as an MxModel object.
}

\usage{
emxRegressionModel(model, data, type='Steven', run, ...)
emxModelRegression(model, data, type='Steven', run, ...)
}


\arguments{
   \item{model}{formula.  See Details.}
   \item{data}{data used for the model}
   \item{run}{logical. Whether to run the model before returning.}
   \item{type}{character.  Either 'Steven' or 'Joshua'.  See Details.}
   \item{...}{Further named arguments to be passed to \code{lm} for the formula}
}

\details{
The \code{model} argument is a formula identical to what is used in \code{\link{lm}}.

The \code{type} argument switches the kind of regression model that is specified.  When there are no missing data, the two versions will estimate the same regression parameters but \code{type='Steven'} will estimate addition parameters that are not estimated by \code{type='Joshua'}.  The \code{type='Steven'} model is due to Steven Boker and many others.  It estimates more parameters than a typical regression analysis and has a different set of assumptions.  More exactly, \code{type='Steven'} models the outcome and all of the predictors as a multivariate Normal distribution.  By contrast, \code{type='Joshua'} is due to Joshua Pritikin and exactly replicates the typical regression model with its usual assumptions.  In particular, \code{type='Joshua'} models the regression residual as a univariate Normal distribution. Predictors are assumed to have no measurement error (see Westfall & Yarkoni, 2016).

The benefit of \code{type='Steven'} is that it handles missing data with
full-information maximum likelihood (FIML; Enders & Bandalos, 2001), at the cost of using a different model with different assumptions from ordinary least squares regression.  The benefit of \code{type='Joshua'} is that it exactly replicates regression as a maximum likelhood model, at the cost of having the same weakness in terms of missing data as OLS regression.

}

\value{
    An MxModel.
}

\seealso{
\link{lm}
}

\references{
Enders, C. K. & Bandalos, D. L. (2001). The relative performance of full
information maximum likelihood estimation for missing data in structural
equation models. <i>Structural Equation Modeling, 8</i>(3), 430-457.

Westfall, J. & Yarkoni, T. (2016). Statistically controlling for
confounding constructs is harder than you think. <i>PLoS ONE,
11</i>(3). doi:10.1371/journal.pone.0152719 }

\examples{
   
# Example
require(EasyMx)
data(myRegDataRaw)
myrdr <- myRegDataRaw
myrdr[1, 4] <- NA

run <- emxRegressionModel(y~1+x*z, data=myrdr, run=TRUE)
summary(run)

summary(lm(y~1+x*z, data=myrdr))



}
