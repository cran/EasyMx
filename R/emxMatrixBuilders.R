#------------------------------------------------------------------------------
matrix2path <- function(x, arrows=1){
	mxPath(
		from=rep(colnames(x@values), each=nrow(x)),
		to=rep(rownames(x@values), times=ncol(x)),
		values=c(x@values),
		free=c(x@free),
		labels=c(x@labels),
		arrows=arrows)
}

emxLoadings <- function(x, values=.8, free=TRUE, path=FALSE){
	ret <- OpenMx::mxMatrix(
		'Full',
		nrow=length(unique(unlist(x))),
		ncol=length(x),
		values=0,
		free=FALSE,
		dimnames=list(unique(unlist(x)), names(x)),
		name='Loadings'
	)
	for(i in 1:length(x)){
		ret@values[x[[names(x)[i]]], names(x)[i]] <- values
		ret@free[x[[names(x)[i]]], names(x)[i]] <- free
		ret@labels[x[[names(x)[i]]], names(x)[i]] <- paste('Load', names(x)[i], x[[names(x)[i]]], sep='')
	}
	if(path) { ret <- matrix2path(ret)}
	return(ret)
}
# Possible alternative interface
#  emxLoadings(x1+x2+x3~F1, x4+x5+x6~F2)

emxResiduals <- function(x, values=.2, free=TRUE, lbound=NA, ubound=NA, path=FALSE, type='unique'){
	if(type=='unique'){
		lab <- paste('Resid', x, sep='')
	} else if(type=='identical'){
		lab <- 'Resid'
	}
	ret <- OpenMx::mxMatrix('Diag', length(x), length(x), free, values, labels=lab, lbound=lbound, ubound=ubound, dimnames=list(x, x), name='Residuals')
	if(path) { ret <- matrix2path(ret)}
	return(ret)
}

emxCovariances <- function(x, values, free, path=FALSE, type, name='Variances'){
	theTypes <- c('independent', 'full', 'corr')
	if(!(type  %in% theTypes)){
		stop(paste("'type' argument not valid.  Valid types are", omxQuotes(theTypes)))
	}
	if(type=='independent'){
		return(emxResiduals(x, values, free, path))
	} else if(type=='full'){
		phval <- matrix(.5, nrow=length(x), ncol=length(x))
		diag(phval) <- 1
		phval <- phval[lower.tri(phval, diag=TRUE)]
		phfre <- TRUE
		phlab <- outer(x, x, paste, sep='')
		diag(phlab) <- paste('Var', diag(phlab), sep='')
		phlab[lower.tri(phlab)] <- paste('Cov', phlab[lower.tri(phlab)], sep='')
		phlab <- phlab[lower.tri(phlab, diag=TRUE)]
		ret <- OpenMx::mxMatrix('Symm', length(x), length(x), phfre, phval, labels=phlab, dimnames=list(x, x), name=name)
	} else if(type=='corr'){
		phval <- matrix(.5, nrow=length(x), ncol=length(x))
		diag(phval) <- 1
		phval <- phval[lower.tri(phval, diag=TRUE)]
		phfre <- phval!=1
		phlab <- outer(x, x, paste, sep='')
		diag(phlab) <- paste('Var', diag(phlab), sep='')
		phlab[lower.tri(phlab)] <- paste('Cov', phlab[lower.tri(phlab)], sep='')
		phlab <- phlab[lower.tri(phlab, diag=TRUE)]
		ret <- OpenMx::mxMatrix('Symm', length(x), length(x), phfre, phval, labels=phlab, dimnames=list(x, x), name=name)
	}
	if(path) { ret <- matrix2path(ret)}
	return(ret)
}


emxMeans <- function(x, values=0, free=TRUE, path=FALSE, type='saturated', name, column=TRUE, labels){
	if(missing(name)){name <- 'Means'}
	if(type=='saturated'){
		lab <- paste('Mean', x, sep='')
	} else if(type=='equal'){
		lab <- paste('M', sep='')
	} else if(type=='twin'){
		lab <- paste('M', x[1:(length(x)/2)], sep='')
	} else if(type=='special'){
		#TODO add error check
		lab <- labels
	}
	if(column) ret <- OpenMx::mxMatrix('Full', length(x), 1, free, values, labels=lab, dimnames=list(x, NULL), name=name)
	else ret <- OpenMx::mxMatrix('Full', 1, length(x), free, values, labels=lab, dimnames=list(NULL, x), name=name)
	if(path) { ret <- matrix2path(ret)}
	return(ret)
}

is.binary <- function(x){
	is.ordered(x) && length(levels(x)) == 2
}

# ordinalCols may be a 1. logical of the same length as the number of columns in data
#  2. character indicating the names of the variables in data
#  3. numeric indicating which columns in the data are ordinal
# TODO: Why not detect which columns are ordinal using is.factor?
emxThresholds <- function(data, ordinalCols=rep(TRUE, ncol(data))){
	if(length(ordinalCols) <= 0){
		stop('You have not specified any ordinal columns.')
	}
	if(is.character(ordinalCols)){
		if(!all( ordinalCols %in% names(data))){
			stop('Some of the ordinal columns requested on not variables in the data.')
		}
		ordinalCols <- match(ordinalCols, names(data))
	}
	if(is.numeric(ordinalCols)){
		tmp <- rep(FALSE, ncol(data))
		tmp[ordinalCols] <- TRUE
		ordinalCols <- tmp
	}
	if(length(ordinalCols) != ncol(data)){
		stop('Weirdness.  You have a difference number of ordinal columns and data columns.')
	}
	numVar <- ncol(data)
	varnam <- names(data)
	ordnam <- colnames(data)[ordinalCols]
	ordinalLevels <- lapply(data[,ordinalCols], levels)
	numOrdinal <- sum(ordinalCols)
	numOrdinalLevels <- sapply(ordinalLevels, length)
	isBinary <- numOrdinalLevels %in% 2
	binnam <- ordnam[isBinary]
	numBinary <- sum(isBinary)
	maxLevels <- max(numOrdinalLevels)
	numThresholds <- maxLevels-1
	thrdnam <- paste(rep(ordnam, each=numThresholds), 'ThrDev', 1:numThresholds, sep='')
	unitLower <- OpenMx::mxMatrix("Lower", numThresholds, numThresholds, values=1, free=FALSE, name="unitLower")
	thrfre <- matrix(NA, numThresholds, numOrdinal)
	thrval <- matrix(NA, numThresholds, numOrdinal)
	for(i in 1:numOrdinal){
		#if(isBinary[i]){
		#	thrfre[,i] <- rep(FALSE, numThresholds)
		#	thrval[,i] <- c(0, rep(.2, numThresholds-1))
		#} else {
			thrfre[,i] <- c(rep(TRUE, numOrdinalLevels[i]-1), rep(FALSE, numThresholds - (numOrdinalLevels[i]-1)))
			thrval[,i] <- rep(.2, numThresholds)
		#}
	}
	thresholdDeviations <- OpenMx::mxMatrix("Full", 
			name="thresholdDeviations", nrow=numThresholds, ncol=numOrdinal,
			values=thrval,
			free = thrfre,
			labels=thrdnam,
			lbound = rep( c(-Inf,rep(.01, (numThresholds-1))) , numOrdinal), # TODO adjust increment value
			dimnames = list(c(), varnam[ordinalCols]),
					)
	saturatedThresholds <- OpenMx::mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix", dimnames=dimnames(thresholdDeviations))
	ret <- list(unitLower, thresholdDeviations, saturatedThresholds)
	if(any(isBinary)){
		Iblock <- diag(1, numBinary)
		colnames(Iblock) <- binnam
		binaryFilterValues <- Iblock
		if (numVar-numBinary > 0) {
			Zblock <- matrix(0, nrow=numBinary, ncol=numVar-numBinary)
			colnames(Zblock) <- varnam[!(varnam %in% binnam)]
			binaryFilterValues <- cbind(binaryFilterValues, Zblock)
		}
		binaryFilterValues <- binaryFilterValues[,varnam]
		BinaryVarianceFilteringMatrix <- NULL  # avoid CRAN check warning
		binaryFilter <- OpenMx::mxMatrix('Full', nrow=numBinary, ncol=numVar, values=binaryFilterValues, free=FALSE, name='BinaryVarianceFilteringMatrix')
		BinaryVarianceFilteringAlgebra <- NULL  # avoid CRAN check warning
		satCov <- NULL
		binaryAlgebraSat <- OpenMx::mxAlgebra(
			BinaryVarianceFilteringMatrix %*% diag2vec(satCov), name='BinaryVarianceFilteringAlgebra')
		BinaryConstantVectorOfOnes <- NULL  # avoid CRAN check warning
		binaryConstant <- OpenMx::mxMatrix('Full', nrow=numBinary, ncol=1, values=1, free=FALSE, name='BinaryConstantVectorOfOnes')
		binaryConstraint <- OpenMx::mxConstraint(
			BinaryConstantVectorOfOnes == BinaryVarianceFilteringAlgebra, name='BinaryVarianceConstraint')
		#ret <- c(ret, list(binaryFilter, binaryAlgebraSat, binaryConstant, binaryConstraint))
	}
	return(ret)
}


#------------------------------------------------------------------------------
