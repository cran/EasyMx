#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2015-06-01
# Filename: emxExamples.R
# Purpose: Show some examples of the emx* functions and user interface
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Source the definitions
require(OpenMx)

require(emx)


#------------------------------------------------------------------------------
# Factor Model Examples

# Example
data(myFADataRaw)
xmap <- list(F1=paste0('x', 1:6), F2=paste0('y', 1:3), F3=paste0('z', 1:3))
mod <- emxFactorModel(xmap, data=myFADataRaw, run=TRUE)

# Example
data(jointdata)
xmap <- list(F=names(jointdata))
mod <- emxFactorModel(xmap, data=jointdata, run=TRUE, ordinal=paste0('z', c(2, 4, 5)))
sat <- mxRefModels(mod, run=TRUE)
summary(mod, refModels=sat)


#Example
require(lavaan)

HS.model <- ' visual =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed =~ x7 + x8 + x9 '
fit <- cfa(HS.model, data=HolzingerSwineford1939)
summary(fit, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)
#inspect(fit, "free")
#inspect(fit, "start")
#inspect(fit, "rsquare")
#inspect(fit, "fit")
#fitted.values(fit)
#coef(fit)
#resid(fit, type="normalized")
modindices(fit)

emap <- list(visual=paste0('x', 1:3), textual=paste0('x', 4:6), speed=paste0('x', 7:9))
efit <- emxFactorModel(emap, data=HolzingerSwineford1939, run=TRUE)



# Example
# OpenMx Frontpage
data(demoOneFactor)
loads <- list(G=names(demoOneFactor))
fit <- emxFactorModel(model=loads, data=demoOneFactor, name='One Factor', run=TRUE)
summary(fit)


#------------------------------------------------------------------------------
# Growth Model Examples

# Example
data(myLongitudinalData)
mod <- emxGrowthModel('Linear', data=myLongitudinalData, use=names(myLongitudinalData), run=TRUE)
c1 <- coef(mod)
# can also do 'Intercept', 'Quadratic', etc., and numbers 0:100+

myLongitudinalData$t0 <- 0
myLongitudinalData$t1 <- 1
myLongitudinalData$t2 <- 2
myLongitudinalData$t3 <- 3
myLongitudinalData$t4 <- 4

mod <- emxGrowthModel('Linear', data=myLongitudinalData, use=names(myLongitudinalData)[1:5], run=TRUE, times=paste0('t', 0:4))

omxCheckCloseEnough(c1, coef(mod))

#------------------------------------------------------------------------------
# Regression with FIML

# Example


data(myRegDataRaw)
myrdr <- myRegDataRaw
myrdr[1, 4] <- NA

run <- emxRegressionModel(y~1+x*z, data=myrdr, run=TRUE)
summary(run)

summary(lm(y~1+x*z, data=myrdr))




#------------------------------------------------------------------------------
# Cholesky model version 1
require(OpenMx)

#Import Data

data(twinData)
twinVar = names(twinData)
selVars <- c('ht1', 'bmi1','ht2','bmi2');  # pick out variables to be modeled (in this case two), for twin 1 and then twin 2
mzdzData <- subset(twinData, zyg %in% c(1, 3), c(selVars, 'zyg'))
# assumes MZ F pairs are coded 1
# assumes DZ f pairs are coded 3
mzdzData$RCoef <- c(1, NA, .5)[mzdzData$zyg]


nVar = length(selVars)/2; # number of dependent variables ** per

x <- paste0('x', 1:nVar)
amat <- emxCholeskyVariance(x, 'A')
cmat <- emxCholeskyVariance(x, 'C')
emat <- emxCholeskyVariance(x, 'E')
ahmat <- emxRelatednessMatrix(2, c(1, .5, 1), labels=c(NA, 'data.RCoef', NA), name='AH')
chmat <- emxRelatednessMatrix(2, c(1, 1, 1), name='CH')
ehmat <- emxRelatednessMatrix(2, c(1, 0, 1), name='EH')
ab <- emxKroneckerVariance(ahmat, amat[[2]], 'AB')
cb <- emxKroneckerVariance(chmat, cmat[[2]], 'CB')
eb <- emxKroneckerVariance(ehmat, emat[[2]], 'EB')
totalVar <- mxAlgebra(AB + CB + EB, 'V', dimnames=list(selVars, selVars))
totalMean <- emxMeans(selVars, type='twin')
expect <- mxExpectationNormal(totalVar$name, totalMean$name)
fitfun <- mxFitFunctionML()

comlist <- list(amat, cmat, emat, ahmat, chmat, ehmat, ab, cb, eb, totalVar, totalMean, expect, fitfun)

model <- mxModel('model', comlist, mxData(mzdzData, 'raw'))
run1 <- mxRun(model)



#------------------------------------------------------------------------------
# Cholesky model version 2 (higher level)

require(OpenMx)
data(twinData)
twinVar = names(twinData)
selVars <- c('ht1', 'bmi1','ht2','bmi2')
mzdzData <- subset(twinData, zyg %in% c(1, 3), c(selVars, 'zyg'))
mzdzData$RCoef <- c(1, NA, .5)[mzdzData$zyg]
nVar = length(selVars)/2
x <- paste0('x', 1:nVar)

acomp <- emxCholeskyComponent(x, 'A', hvalues=c(1, .5, 1), hlabels=c(NA, 'data.RCoef', NA))
ccomp <- emxCholeskyComponent(x, 'C', hvalues=c(1, 1, 1))
ecomp <- emxCholeskyComponent(x, 'E', hvalues=c(1, 0, 1))
totalVar <- mxAlgebra(AKron + CKron + EKron, 'V', dimnames=list(selVars, selVars))
totalMean <- emxMeans(selVars, type='twin')
expect <- mxExpectationNormal(totalVar$name, totalMean$name)
fitfun <- mxFitFunctionML()

comlist <- c(acomp, ccomp, ecomp, list(totalVar, totalMean, expect, fitfun))

model <- mxModel('model', comlist, mxData(mzdzData, 'raw'))
run2 <- mxRun(model)


#------------------------------------------------------------------------------
# Cholesky model version 3 (highest level)

require(OpenMx)
data(twinData)
twinVar = names(twinData)
selVars <- c('ht1', 'bmi1','ht2','bmi2')
mzdzData <- subset(twinData, zyg %in% c(1, 3), c(selVars, 'zyg'))
mzdzData$RCoef <- c(1, NA, .5)[mzdzData$zyg]

run3 <- emxTwinModel(model='Cholesky', relatedness='RCoef',
	data=mzdzData, use=selVars, run=TRUE, name='TwCh')


#------------------------------------------------------------------------------
# Genetic factor model

data(twinData)
twinVar = names(twinData)
selVars <- c('ht1', 'bmi1','wt1','htwt1','ht2','bmi2','wt2','htwt2');  # pick out variables to be modeled (in this case two), for twin 1 and then twin 2
mzdzData <- subset(twinData, zyg %in% c(1, 3), c(selVars, 'zyg'))
# assumes MZ F pairs are coded 1
# assumes DZ f pairs are coded 3
mzdzData$RCoef <- c(1, NA, .5)[mzdzData$zyg]
mzdzData[,selVars] <- scale(mzdzData[,selVars])

nVar = length(selVars)/2; # number of dependent variables ** per
x <- selVars[1:4]

acomp <- emxGeneticFactorComponent(x, 'A', hvalues=c(1, .5, 1), hlabels=c(NA, 'data.RCoef', NA))
ccomp <- emxGeneticFactorComponent(x, 'C', 1, FALSE, hvalues=c(1, 1, 1))
ecomp <- emxGeneticFactorComponent(x, 'E', .5, TRUE, hvalues=c(1, 0, 1))
rcomp <- emxResiduals(selVars, lbound=1e-7)
totalVar <- mxAlgebra(AKron + CKron + EKron + Residuals, 'V', dimnames=list(selVars, selVars))
totalMean <- emxMeans(selVars, type='twin', column=FALSE)
expect <- mxExpectationNormal(totalVar$name, totalMean$name)
fitfun <- mxFitFunctionML()

comlist <- c(acomp, ccomp, ecomp, list(rcomp, totalVar, totalMean, expect, fitfun))

model <- mxModel('model', comlist, mxData(mzdzData, 'raw'))
run <- mxRun(model)





