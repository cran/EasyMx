

emxMixtureModel <- function(model, data, run=FALSE, p=NA, ...){
	models <- model
	# TODO handle data better
	data <- mxData(data, 'raw')
	modelNames <- sapply(models, slot, 'name')
	# Check that all models have ML fit functions
	for(i in 1:length(models)){
		if(!(class(models[[i]]$fitfunction) %in% 'MxFitFunctionML')){
			stop('Cannot proceed because model ', OpenMx::omxQuotes(modelNames[i]), ' does not have an ML fit function.')
		}
		models[[i]]$fitfunction <- OpenMx::mxFitFunctionML(vector=TRUE, rowDiagnostics=models[[i]]$fitfunction$rowDiagnostics)
		models[[i]]$data <- NULL
	}
	model <- OpenMx::mxModel(model='MixtureModel', models, data)
	if(length(p)==1 && is.na(p)){
		p1 <- OpenMx::mxMatrix('Full', nrow=length(models), ncol=1, values=1, lbound=1e-6, free=c(FALSE, rep(TRUE, length(models)-1)), name='theSmallPMatrix', labels=paste0('unscaledProportionForModel', modelNames))
		theSmallPMatrix <- NULL
		p2 <- OpenMx::mxAlgebra((1/sum(theSmallPMatrix)) %x% theSmallPMatrix, name='theScaledPMatrix')
		p <- 'theScaledPMatrix'
		model <- OpenMx::mxModel(model, p1, p2, ...)
	} else {
		model <- OpenMx::mxModel(model, ...)
	}
	theAlg1 <- OpenMx::mxAlgebraFromString(paste('cbind(', paste(modelNames, '.fitfunction', sep='', collapse=', '), ') %*% ', p, sep=''), name='theMixtureFitVector')
	theMixtureFitVector <- NULL
	theAlg2 <- OpenMx::mxAlgebra(-2*sum(log(theMixtureFitVector)), name='theMixtureFit')
	theFit <- OpenMx::mxFitFunctionAlgebra('theMixtureFit')
	model <- OpenMx::mxModel(model, theAlg1, theAlg2, theFit)
	if(run){
		model <- mxRun(model)
	}
	return(model)
}

emxModelMixture <- emxMixtureModel

#v2 <- OpenMx::mxMixtureModel(list(class1, class2), data=dataRaw)
#v2r <- OpenMx::mxRun(v2)
#summary(v2r)

