#------------------------------------------------------------------------------
# Variance components

emxVarianceComponents <- function(model, data, run){
	
	if(run==TRUE){
		model <- OpenMx::mxRun(model)
	}
	return(model)
}


#model=list(variance=, relatedness=)

#--------------------------------------
# Mid-level Functions

emxCholeskyVariance <- function(x, name, values=.8, free=TRUE){
	# TODO process the actual names of the variables given in x
	# These can be the dimnames of the OpenMx::mxMatrix and OpenMx::mxAlgebra
	nvar <- length(x)
	sqrtName <- paste0('sqrt', name)
	labs <- paste0(sqrtName, outer(1:nvar, 1:nvar, paste0)[lower.tri(matrix(0, nvar, nvar), diag=TRUE)])
	lowerm <- OpenMx::mxMatrix('Lower', nvar, nvar, free, values, labels=labs, name=sqrtName)
	diag(lowerm$lbound) <- 1e-6
	algText <- paste0('mxAlgebra(', sqrtName,' %*% t(', sqrtName, '), name="', name, '")')
	fullm <- eval(parse(text=algText))
	return(list(lowerm, fullm))
}

emxRelatednessMatrix <- function(nvar, values, labels, name='h'){
	if(missing(labels)){labels <- paste0(name, 1:(nvar*(nvar+1)/2))}
	mxMatrix('Symm', nvar, nvar, free=FALSE, values=values, labels=labels, name=name)
}

emxKroneckerVariance <- function(h, v, name){
	if(missing(name)){name <- paste0('kron', v$name)}
	algText <- paste0('mxAlgebra(', h$name, ' %x% ', v$name, ', name=name)')
	ret <- eval(parse(text=algText))
	return(ret)
}

emxGeneticFactorVariance <- function(x, name, values=.8, free=TRUE, lbound=NA, ubound=NA){
	nvar <- length(x)
	sqrtName <- paste0('inner', name)
	labs <- paste0(sqrtName, 1:nvar)
	lowerm <- OpenMx::mxMatrix('Full', nvar, 1, free, values, labels=labs, name=sqrtName, lbound=lbound, ubound=ubound)
	algText <- paste0('mxAlgebra(', sqrtName,' %*% t(', sqrtName, '), name=name)')
	fullm <- eval(parse(text=algText))
	return(list(lowerm, fullm))
}



#--------------------------------------
# Slightly higher-level Functions


emxCholeskyComponent <- function(x, xname, h=2, hname=paste0('H', xname), hvalues, hlabels){
	xmat <- emxCholeskyVariance(x, xname)
	hmat <- emxRelatednessMatrix(h, values=hvalues, labels=hlabels, name=hname)
	kmat <- emxKroneckerVariance(hmat, xmat[[2]], paste0(xname, 'Kron'))
	return(list(xmat[[1]], xmat[[2]], hmat, kmat))
}

emxGeneticFactorComponent <- function(x, xname, xvalues=.8, xfree=TRUE, xlbound=NA, xubound=NA, h=2, hname=paste0('H', xname), hvalues, hlabels){
	xmat <- emxGeneticFactorVariance(x, xname, xvalues, xfree, xlbound, xubound)
	hmat <- emxRelatednessMatrix(h, values=hvalues, labels=hlabels, name=hname)
	kmat <- emxKroneckerVariance(hmat, xmat[[2]], paste0(xname, 'Kron'))
	return(list(xmat[[1]], xmat[[2]], hmat, kmat))
}

#TODO
# Genetic Factor model
# A  = la %*% t(la)
# C  = lc %*% t(lc)
# E  = le %*% t(le)
# Cov = A+C+E
#       h*A+C  A+C+E
#
#  la, lc, le are Nx1
#
#
# Add Independent pathway model component function
# It's basically a genetic factor model for each A/C/E component
#  plus a uniquenesses for each component.
# A  = la %*% t(la) + ua
# C  = lc %*% t(lc) + uc
# E  = le %*% t(le) + ue
# Cov = A+C+E
#       h*A+C  A+C+E
#
#  la, lc, le are Nx1
#  ua, uc, ue are NxN diagonal
#
#
# Common Pathway model component function
# A  = lf %*% la %*% t(la) %*% lf + ua
# C  = lf %*% lc %*% t(lc) %*% lf + uc
# E  = lf %*% le %*% t(le) %*% lf + ue
# Cov = A+C+E
#       h*A+C  A+C+E
#
#  lf is Nx1
#  la, lc, le are 1x1
#  ua, uc, ue are NxN diagonal
#
#



#--------------------------------------
# Very-high level Functions

#model=Factor, Cholesky, IndependentPathway, CommonPathway
emxTwinModel <- function(model, relatedness, data, run=FALSE, use, name='model'){
	#TODO Add processing for the model argument
	#TODO How to specify subsets of A, D, C, E components
	#TODO Parse relatedness argument
		#When string, interpret as name of variables in data set giving coeficient of relatedness
		#When single number, interpret as numeric value of the coefficient of relatedness
	rlen <- length(relatedness)
	if(rlen == 1){
		if(is.numeric(relatedness)){
			hval <- c(1, relatedness, 1)
			hlab <- c(NA, NA, NA)
		} else if(is.character(relatedness)){
			hval <- c(1, .5, .1)
			hlab <- c(NA, paste0('data.', relatedness), NA)
		}
	} else if(rlen == nrow(data)){
		data <- cbind(data, RelatednessCoefficient=relatedness)
		hval <- c(1, .5, .1)
		hlab <- c(NA, paste0('data.', 'RelatednessCoefficient'), NA)
	}
	x <- paste0('x', 1:(length(use)/2))
	acomp <- emxCholeskyComponent(x, 'A', hvalues=hval, hlabels=hlab)
	ccomp <- emxCholeskyComponent(x, 'C', hvalues=c(1, 1, 1))
	ecomp <- emxCholeskyComponent(x, 'E', hvalues=c(1, 0, 1))
	AKron <- NULL
	CKron <- NULL
	EKron <- NULL
	totalVar <- OpenMx::mxAlgebra(AKron + CKron + EKron, 'V', dimnames=list(use, use))
	totalMean <- emxMeans(use, type='twin')
	expect <- OpenMx::mxExpectationNormal(totalVar$name, totalMean$name)
	fitfun <- OpenMx::mxFitFunctionML()
	
	comlist <- c(acomp, ccomp, ecomp, list(totalVar, totalMean, expect, fitfun))
	
	model <- OpenMx::mxModel(model=name, comlist, OpenMx::mxData(data, 'raw'))
	if(run){
		model <- OpenMx::mxRun(model)
	}
	return(model)
}

#emxCholeskyModel(model=c('A','C', 'E'), relatedness=c('standard', 'RCoef'), data=mzdzData)

emxModelTwin <- emxTwinModel

