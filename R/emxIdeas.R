##------------------------------------------------------------------------------
#require(OpenMx)
#data <- read.csv("data.csv")
#result <- emxFactorModel(model=list(F=names(data)),
#	data=data, name='A Great Name', run=TRUE)

#emxFactorModel(model, data, name, run, identification, use, ordinal)
## model is named list
##list(F=names(data) creates a one-factor model with F as the factor and names(data) as the indicators
##
##data is a data.frame or matrix
##run is logical, whether or not to run the model once created
##identification is a character, how to identifiy the model (fix factor loading, fix mean/variance of factor)
##use is a character vector, the names of the variable to use.
##  defaults to the variables used in the model
##ordinal is a character vector, the names of the ordinal variables
##
## maybe have LID (latent identification) and OID (ordinal identification)


##------------------------------------------------------------------------------
#require(OpenMx)
#data <- read.csv("data.csv")
#result <- emxGrowthModel(model="Linear", data=data)

#emxGrowthModel(model, data, name, run, identification, use, ordinal)
## model is a single character string
## "Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic", 0:100,
## "Latent" (latent basis)


##------------------------------------------------------------------------------
#require(OpenMx)
#data <- read.csv("data.csv")
#result <- emxLDEModel(model=)

#emxLDEModel(model, data, name, run, identification, use, ordinal)
## need to give the window size (embedding dimension), the order of polynomial approximation,
## the order of the derivative


##------------------------------------------------------------------------------
#require(OpenMx)
#data <- read.csv("data.csv")
#result <- emxLDEModel(model=)


##------------------------------------------------------------------------------
#require(OpenMx)
#data <- read.csv("data.csv")
#result <- emxACEModel(model=)

## Martin & Eaves (1977)
#result <- emxGeneticFactorModel(model, data, name, run=FALSE, identification, use, ordinal)

## Cholesky
#result <- emxCholeskyModel(model=)

## independent pathways
#result <- emxIPModel(model=)

## common pathways
#result <- emxCPModel(model=)

