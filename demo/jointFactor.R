#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-09-08
# Filename: jointFactorBug.R
# Purpose: Originaly was trying to fit a WLS model and then get factor
#  scores, but then found LISREL WLS seems to fail.  That's the real
#  problem.

#------------------------------------------------------------------------------
# Load packages and data

require(OpenMx)
require(EasyMx)


data(jointdata)
# specify ordinal columns as ordered factors
jointdata[,c(2,4,5)] <- mxFactor(jointdata[,c(2,4,5)], 
	levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))


#------------------------------------------------------------------------------
# ML version of model
mlMod <- emxFactorModel(list(F=names(jointdata)), jointdata)
mlMod$Loadings$lbound <- 0
mlFit <- mxRun(mlMod)
# takes a few seconds, but works

mxEval(thresholdMatrix, mlFit)
mxEval(thresholdDeviations, mlFit)


#------------------------------------------------------------------------------
# WLS version of model

wlsMod <- mxModel(mlMod, mxDataWLS(jointdata), mxFitFunctionWLS())
wlsFit <- mxRun(wlsMod)


#------------------------------------------------------------------------------
# Compare estimates

cbind(ML=coef(mlFit), WLS=coef(wlsFit))



#------------------------------------------------------------------------------
# End
