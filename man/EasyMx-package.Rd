\name{EasyMx-package}
\alias{EasyMx-package}
\alias{EasyMx}
\docType{package}

\title{EasyMx: Easy modeling in OpenMx}

\description{EasyMx is a package for extended structural equation modeling. It is built as a higher-level frontend on top of OpenMx.  It is intended as an Easy introduction to OpenMx: Easy Mx.  Try the example below.}

\details{
All of the functions in the EasyMx package create OpenMx objects. These are most often MxMatrix, MxAlgebra, or MxModel objects.  The primary difference between EasyMx and OpenMx is design philosophy.  OpenMx has its foundation in WYSIWID: What you say is what it does.  This requires the user to be very explicit.  The EasyMx package is not as strong or flexible as OpenMx, but it places less burden on the user.  Many decisions are made automatically for the user.  Some of them are modifiable within EasyMx; for others the user is encouraged to use OpenMx, where nearly everything is modifiable.

The package is broadly divided into two styles of functions: matrix builders and model builders.

The matrix builder functions are utilities for building common structural equation model matrices.  These include \code{\link{emxLoadings}} for factor loadings, \code{\link{emxResiduals}} for residual variances, \code{\link{emxCovariances}} for latent or manifest covariances, \code{\link{emxMeans}} for means and intercepts matrices, and \code{\link{emxThresholds}} for thresholds matrices when ordinal data are involved.

The model builder functions are higher-level utilities for building common kinds of structural equation models.  The model builders often call several matrix builders.  The model builders include \code{\link{emxFactorModel}} for (multiple) factor models, \code{\link{emxGrowthModel}} for latent growth curve models, and \code{\link{emxRegressionModel}} for full-information likelihood estimation of regression for observed variables.

A third category of functions encompasses special functions for behavior genetics modeling.  Some of these functions are matrix builders, and others are model builders.  The lowest-level functions for behavior genetics are \code{\link{emxCholeskyVariance}}, \code{\link{emxGeneticFactorVariance}}, \code{\link{emxRelatednessMatrix}}, and \code{\link{emxKroneckerVariance}}.

A higher-level set of behavior genetics matrix builders create all the matrices and algebraic statements needed for e.g. the A component of an ACE model.  These functions are \code{\link{emxCholeskyComponent}} and \code{\link{emxGeneticFactorComponent}}.

The highest-level of behavior genetics functions builds some basic twin models.  The primary function for this is \code{\link{emxTwinModel}}.

Finally, a mixture model helper is provided: \code{\link{emxMixtureModel}}.
}


\examples{

# Make and run a one factor model
\dontrun{
require(EasyMx)

data(demoOneFactor)
fmod <- list(G=names(demoOneFactor))
fit1 <- emxFactorModel(fmod, demoOneFactor, run=TRUE)
summary(fit1)
}
}

