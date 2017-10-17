




#------------------------------------------------------------------------------
buildItemFactorModel <- function(model, name, bdata, use) {
	if (!requireNamespace("rpf", quietly=TRUE)) stop("Please install package 'rpf' to use item factor models")
	spec <- lapply(use, function(x) rpf::rpf.grm(outcomes=length(levels(bdata$observed[,x]))))
	names(spec) <- use
	
	imat <- mxMatrix(name="item", values=mxSimplify2Array(lapply(spec, rpf::rpf.rparam)))
	imat$free[,] <- !is.na(imat$values)
	rownames(imat)[1:length(model)] <- names(model)
	for (fx in 1:length(model)) {
		noLoading <- setdiff(use, model[[fx]])
		imat$values[fx,noLoading] <- 0
		imat$free[fx,noLoading] <- FALSE
	}
	
	qpts <- 49L
	if (length(model) == 2) qpts <- 31L
	if (length(model) == 3) qpts <- 21L
	
	plan <- mxComputeSequence(list(mxComputeEM(
		'expectation', 'scores',
		mxComputeNewtonRaphson(verbose=0L), information="oakes1999",
	infoArgs=list(fitfunction='fitfunction')),
	mxComputeOnce('fitfunction', 'gradient'),
	mxComputeHessianQuality(),
	mxComputeStandardError(),
	mxComputeReportDeriv()))
	
	mxModel(name=name, imat, bdata,
			mxExpectationBA81(ItemSpec=spec, qpoints=qpts),
			mxFitFunctionML(),
			plan)
}

buildLisrelFactorModel <- function(model, name, bdata, data, use, ordinalCols) {
	latents <- names(model)
	mmat <- emxMeans(x=use, free=!ordinalCols)
	lmat <- emxLoadings(x=model)
	rmat <- emxResiduals(x=use, free=!ordinalCols, values=1)
	ka <- emxMeans(x=latents, free=FALSE, name='LatentMeans')
	ph <- emxCovariances(x=latents, type='corr', name='LatentVariances')
	bmodel <- OpenMx::mxModel(name=name,
		lmat, rmat, mmat, ka, ph,
		bdata,
		mxExpectationLISREL(
			LX=slot(lmat, 'name'),
			TD=slot(rmat, 'name'),
			TX=slot(mmat, 'name'),
			KA=slot(ka, 'name'),
			PH=slot(ph, 'name')),
		mxFitFunctionML()
	)
	if(any(ordinalCols)){
		threshList <- emxThresholds(data, ordinalCols)
		bmodel <- OpenMx::mxModel(bmodel, threshList,
			mxExpectationLISREL(
			LX=slot(lmat, 'name'),
			TD=slot(rmat, 'name'),
			TX=slot(mmat, 'name'),
			KA=slot(ka, 'name'),
			PH=slot(ph, 'name'),
			thresholds=slot(threshList[[3]], 'name'))
		)
	}
	bmodel
}

emxFactorModel <- function(model, data, name, run=FALSE, identification, use, ordinal, ...,
							parameterization=c("lisrel", "ifa"), weight = as.character(NA)) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("emxFactorModel does not accept values for the '...' argument")
	}
	if(missing(name)){name <- 'Model'}
	if(missing(use)){
		use <- unlist(model)
		use <- use[!duplicated(use)]
	}
	numVar <- length(use)
	parameterization <- match.arg(parameterization)
	# TODO Write more general data processing module
	if( nrow(data) == ncol(data) && all(data == t(data)) ){
		if (parameterization == 'ifa') stop("Only raw data can be used with parameterization='ifa'")
		if (!is.na(weight)) stop("Cannot use row weights with covariance data")
		data <- data[use, use]
		bdata <- OpenMx::mxData(data, 'cov')
	} else {
		weightCol <- NULL
		if (!is.na(weight)) weightCol <- data[[weight]]
		missingCol <- !(use %in% colnames(data))
		if (any(missingCol)) {
			stop(paste("Cannot find column", omxQuotes(use[missingCol]), "in the data"))
		}
		data <- data[,use,drop=FALSE]
		if(missing(ordinal)){
			if (is.data.frame(data)) {
				ordinalCols <- sapply(data, is.ordered)
			} else {
				ordinalCols <- rep(FALSE, numVar)
			}
		} else {
			if(!all(ordinal %in% use)){stop('Some of the ordinal variables are not among those being used.  Either include the ordinal variables in use, or exclude the ordinal variables not being used.')}
			ordinalCols <- rep(FALSE, ncol(data))
			ordinalCols[match(ordinal, use)] <- TRUE
		}
		if(!any(is.na(data)) && !any(ordinalCols) && is.na(weight)) {
			bdata <- OpenMx::mxData(cov(data), 'cov', means=colMeans(data), numObs=nrow(data))
		} else {
			unconverted <- ordinalCols & !sapply(data, is.ordered)
			if(any(unconverted)) {
				for (cx in which(unconverted)) {
					data[[cx]] <- mxFactor(data[[cx]], levels=sort(unique(data[[cx]])))
				}
			}
			weightedData <- data
			if (!is.na(weight)) weightedData[[weight]] <- weightCol
			bdata <- OpenMx::mxData(weightedData, 'raw', weight=weight)
		}
	}
	if (parameterization == 'lisrel') {
		bmodel <- buildLisrelFactorModel(model, name, bdata, data, use, ordinalCols)
	} else if (parameterization == 'ifa') {
		if (!all(ordinalCols)) stop(paste("All columns must be ordinal for parameterization=",
			omxQuotes(parameterization)))
		if (length(model) > 3L) stop("More than 3 latent factors are not supported for IFA models")
		bmodel <- buildItemFactorModel(model, name, bdata, use)
	} else { stop(parameterization) }
	if(run){return(mxRun(bmodel))}
	return(bmodel)
}

emxModelFactor <- emxFactorModel






#------------------------------------------------------------------------------
# Growth model

GrowthBasisMatrix <- function(numTimes, deltaT=1, order=2, steps=NA){
	if(all(is.na(steps))){
		m <- numTimes
		steps <- (0:(m-1))*deltaT
	}
	L <- outer(steps, 0:order, "^")
	return(L)
}


resolveGrowthModel <- function(x){
	nums <- c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic")
	if(length(x) != 1){
		stop('wtf am I supposed to do with that?')
	}
	if(is.character(x)){
		if( x %in% nums ){
			return( match(x, nums) - 1 )
		} else {stop(
			paste('I do not understand what model is meant by', omxQuotes(x), '. ', 
			'The model argument must be a number or one of the following: ', 
			omxQuotes(nums), sep=''))
		}
	}
	else if( round(x) == x && x >= 0){
		return(x)
	}
	else{stop('I do not understand your model arument.')}
}



emxGrowthModel <- function(model, data, name, run=FALSE, identification, use, ordinal, times){
	if(missing(name)){name <- 'Model'}
	if(missing(use)){stop('You must specify use for growth models')}
	if(missing(times)){times <- 0:(length(use)-1)}
	order <- resolveGrowthModel(model)
	latents <- paste('F', 0:order, sep='')
	manifests <- use
	# TODO Write more general data processing module
	if( nrow(data) == ncol(data) && all(data == t(data)) ){
		bdata <- OpenMx::mxData(data, 'cov')
	} else {
		#if(!any(is.na(data))){
		#	bdata <- OpenMx::mxData(cov(data), 'cov', means=colMeans(data), numObs=nrow(data))
		#} else {
		bdata <- OpenMx::mxData(data, 'raw')
		#}
	}
	mmat <- emxMeans(x=use, values=0, free=FALSE, type='equal')
	if(is.numeric(times)){
		lval <- GrowthBasisMatrix(steps=times, order=order)
		lmat <- OpenMx::mxMatrix('Full', length(use), length(latents), FALSE, values=lval, name='Loadings', dimnames=list(manifests, latents))
		tmat <- NULL
	} else if(is.character(times)){
		tmat <- OpenMx::mxMatrix('Full', length(use), 1, FALSE, values=1, name='Times', labels=paste0('data.', times))
		tstring <- paste0('cbind( ', paste('Times^', 0:order, sep='', collapse=', '), ' )')
		lmat <- OpenMx::mxAlgebraFromString(tstring, name='Loadings', dimnames=list(manifests, latents))
	}
	rmat <- emxResiduals(x=use, type='identical')
	ka <- emxMeans(x=latents, free=TRUE, name='LatentMeans')
	ph <- emxCovariances(x=latents, type='full', name='LatentVariances')
	bmodel <- OpenMx::mxModel(name=name,
		lmat, rmat, mmat, ka, ph, tmat,
		bdata,
		mxExpectationLISREL(
			LX=slot(lmat, 'name'),
			TD=slot(rmat, 'name'),
			TX=slot(mmat, 'name'),
			KA=slot(ka, 'name'),
			PH=slot(ph, 'name')),
		mxFitFunctionML()
	)
	if(run){return(mxRun(bmodel))}
	return(bmodel)
}

emxModelGrowth <- emxGrowthModel




#------------------------------------------------------------------------------
# Regression with FIML

emxRegressionModel <- function(model, data, type='Steven', run, ...){
	theFormu <- model #a formula as passed to lm
	theFrame <- lm(theFormu, data=data, na.action=na.pass, method='model.frame', ...)
	theTerms <- attr(theFrame, 'terms')
	theRespo <- model.response(theFrame, 'numeric')
	theMatri <- model.matrix(theTerms, theFrame)
	mnam <- gsub(":", "_x_", colnames(theMatri))
	colnames(theMatri) <- mnam
	whichInt <- match('(Intercept)', colnames(theMatri))
	if(!is.na(whichInt)){
		theMatri <- theMatri[, -whichInt]
	}
	
	namRespo <- names(theFrame)[1]
	namTerms <- colnames(theMatri)
	namAll <- c(namRespo, namTerms)
	cdat <- cbind(theRespo, theMatri)
	colnames(cdat) <- namAll
	
	theData <- OpenMx::mxData(cdat, 'raw')
	theInter <- OpenMx::mxPath(from='one', to=namRespo, labels='Intercept', values=0, free=!is.na(whichInt))
	residVar <- OpenMx::mxPath(from=namRespo, arrows=2, free=TRUE, values=1, labels=paste('residVar', namRespo, sep='_'))
	
	if (type=='Joshua'){
		namLat <- paste0(namTerms, "Latent")
		structPaths <- OpenMx::mxPath(from='one', to=namLat, free=FALSE, labels=paste0('data.', namTerms))
		regPaths <- OpenMx::mxPath(from=namLat, to=namRespo, labels=paste0(namRespo, '_on_', namTerms))
		model <- OpenMx::mxModel('model', type='RAM', manifestVars=namRespo, latentVars=namLat, theData, theInter, residVar, structPaths, regPaths)
	} else if(type=='Steven'){
		regPaths <- OpenMx::mxPath(from=namTerms, to=namRespo, free=TRUE, values=.8, labels=paste(namRespo, '_on_', namTerms, sep=''))
		predVars <- OpenMx::mxPath(from=namTerms, arrows=2, free=TRUE, values=1, labels=paste('var', namTerms, sep='_'))
		theMeans <- OpenMx::mxPath(from='one', to=namTerms, labels=paste('mean', namTerms, sep='_'))
		model <- OpenMx::mxModel('model', type='RAM', manifestVars=namAll, theInter, regPaths, residVar, theMeans, predVars, theData)
	} else {
		stop(paste('Unknown type argument, ', type, ', in emxRegressionModel', sep=''))
	}
	
	if(run==TRUE){
		model <- OpenMx::mxRun(model)
	}
	return(model)
}


emxModelRegression <- emxRegressionModel




