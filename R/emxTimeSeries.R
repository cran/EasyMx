#------------------------------------------------------------------------------


buildVARModel <- function(use, bdata, name){
	latents <- paste0('Latent_', use)

	mmat <- emxMeans(x=use, free=FALSE, values=0)

	llist <- as.list(use)
	names(llist) <- latents
	lmat <- emxLoadings(x=llist, values=1, free=FALSE)

	rmat <- emxResiduals(x=use, free=FALSE, values=0)

	ka <- emxMeans(x=latents, free=TRUE, name='LatentMeans')

	nlatents <- length(latents)
	ph <- OpenMx::mxMatrix('Diag', nlatents, nlatents, free=TRUE, values=.5, name='DynamicNoiseCov', lbound=1e-6, labels=paste0('DynNoise', latents))
	
	dy <- OpenMx::mxMatrix('Full', nlatents, nlatents,
		values=diag(.2, nrow=nlatents), free=TRUE,
		labels=paste0('ar_', outer(1:nlatents, 1:nlatents, paste0)),
		lbound=-2.5, ubound=2.5,
		name='Dynamics', dimnames=list(latents, latents))
	
	im <- emxMeans(paste0(latents, '_0'), free=FALSE, name='x0')
	ic <- OpenMx::mxMatrix('Diag', nlatents, nlatents, labels=paste0(latents, '_0'), free=FALSE, values=1, name='P0')
	uu <- OpenMx::mxMatrix('Full', 1, 1, values=1, name='u')
	
	bmodel <- OpenMx::mxModel(name=name,
		lmat, rmat, mmat, ka, ph, dy, im, ic, uu,
		bdata,
		mxExpectationStateSpace(
			A=slot(dy, 'name'),
			B=slot(ka, 'name'),
			C=slot(lmat, 'name'),
			D=slot(mmat, 'name'),
			Q=slot(ph, 'name'),
			R=slot(rmat, 'name'),
			x0=slot(im, 'name'),
			P0=slot(ic, 'name'),
			u=slot(uu, 'name')),
		mxFitFunctionML()
	)

}

emxVARModel <- function(model, data, name, run=FALSE, use, ID) {
	model=c('diag', 'full')
	# TODO lags!
	data <- data[,use, drop=FALSE]
	bdata <- OpenMx::mxData(data, 'raw')

	m <- buildVARModel(use=use, bdata=bdata, name=name)
	if(run){return(mxRun(m))} else {return(m)}
}

emxModelVAR <- emxVARModel

# emxLVARModel <- function(model, data, name, run=FALSE, use, ID) {
	
# }

# emxHarveyModel <- function(model, data, name, run=FALSE, use, ID) {
# 	random walk
# 		drift='zero', measurement='fixed'
# 	random walk with drift
# 		drift='constant', measurement='fixed'
# 	local level = random walk with measurement noise
# 		drift='zero', measurement='stochastic'
# 	stochastic trend = random walk with autoregressive drift
# 		drift='AR', measurement='fixed'
# 	local linear trend = random walk with autoregressive drift and measurement noise
# 		drift='AR', measurement='stochastic'
# 	drift = c('zero', 'constant', 'AR')
# 	measurement = c('fixed', 'stochastic')
	
# }

#emxStateSpaceModel <- function(model, data, name, run=FALSE, use, ID, time, ...) {
#	stop('This function is not yet implemented')
#}

emxStateSpaceMixtureModel <- function(model, data, name, run=FALSE, use, ID, time, ...) {
	# model = "VAR", a single MxModel, a list of MxModel objects
	# mix = c('dynamics', 'measurement', 'all')

	if(missing(use)){
		# extract dimnames from all 'C' matrices in model list and use the union of those
	}

	if(is.data.frame(data)){
		# Split data into list by ID
		ids <- unique(data[,ID])
		pdim <- length(ids)
		dsList <- list()
		for(p in 1:pdim){
			dsList[[p]] <- data[data[,ID] == ids[p],]
		}
	} else if(is.list(data)){
		pdim <- length(data)
		ids <- 1:pdim #sapply(sapply(data, '[', ID), '[', i=1)
		dsList <- data
	} else {stop("'data' argument should be a list of data.frame objects or a single data.frame with ID column")}

	# Create a model for each person in each mixture class
	kmods <- if(is.list(model)) model else {stop("'model' argument should be a list of models")}
	kdim <- length(model)
	pkmodels <- list()
	pk <- 1
	for(p in 1:pdim){
	    for(k in 1:kdim){
	        pkmodels[[pk]] <- mxModel(kmods[[k]],
	            name=paste0('Person', ids[p], 'Klass', k),
	            mxData(dsList[[p]], 'raw'),
	            mxFitFunctionML(vector=TRUE))
	        pk <- pk + 1
	    }
	}
	names(pkmodels) <- sapply(pkmodels, slot, 'name')

	# Mixture-specific matrix
	if(kdim == 1) {kprob <- 1} else if(kdim == 2) {kprob <- c(2/3, 1/3)} else {kprob <- c(1/(kdim-1), rep((1 - 1/(kdim-1))/(kdim-1), kdim-1))}
	#c(.5, .25, .25)
	kmat <- mxMatrix('Full', nrow=kdim, ncol=1, values=log(kprob/min(kprob)),
    	labels=paste0('klass_prob', 1:kdim), free=c(rep(TRUE, kdim-1), FALSE), name='K')

	# Create a mixture model for each person that contains the class-person
	#  combined models for that person
	pmodels <- list()
	for(p in 1:pdim){
	    pind <- grep(paste0('Person', ids[p], 'Klass'), names(pkmodels))
	    pmodels[[p]] <- mxModel(model=paste0('Person', ids[p]), pkmodels[pind],
	        kmat,
	        mxExpectationMixture(
	            components=names(pkmodels)[pind],
	            weights='K', scale='softmax'),
	        mxFitFunctionML()
	    )
	}
	names(pmodels) <- sapply(pmodels, slot, 'name')

	# Create a multigroup model composed of all the mixture models for each person
	grpmodel <- mxModel(model='container', pmodels,
	    mxFitFunctionMultigroup(names(pmodels))
	)

	if(run){return(mxRun(grpmodel))} else {return(grpmodel)}
}


emxModelStateSpaceMixture <- emxStateSpaceMixtureModel


emxStateSpaceMixtureClassify <- function(model){
	pdim <- length(model$submodels)
	pknames <- sapply(sapply(model$submodels, slot, name='submodels'), slot, name='name')
	kdim <- length(pknames)/pdim
	kest <- mxEvalByName(paste0(slot(model$submodels[[1]], name='name'), '.K'), model)
	kestp <- exp(kest)/sum(exp(kest))
	# TODO generalize this for other scalings
	# TODO check that model wasRun
	# TODO check that input model is something like a state space mixture model

	# matrix of fit function names for every person-class combination
	pkmat <- matrix(
    	paste0(pknames, '.fitfunction'),
    	nrow=pdim, ncol=kdim, byrow=TRUE)

	# row likelihoods for every person-class-time
	likArray <- apply(pkmat, c(1, 2), mxEvalByName, model=model)
	# dim 1: tdim number of time points
	# dim 2: pdim number of people
	# dim 3: kdim number of classes

	# minus 2 log likelihood for each person in every class
	m2llMat <- apply(likArray, c(2, 3), function(x){ -2*sum(log(x))})
	# dim 1: pdim number of people
	# dim 2: kdim number of classes

	# combine minus 2 log lik with prior prob of each class
	#  in the log scale
	fullM2llMat <- m2llMat + 
	    matrix(-2*log(kestp), nrow=pdim, ncol=kdim, byrow=TRUE)

	# best fitting class for each person
	est_klasses <- apply(fullM2llMat, 1, which.min)
	
	# return object
	ret <- list(estimated_classes=est_klasses,
		joint_m2ll=fullM2llMat, m2ll=m2llMat, likelihood=likArray)
	
	return(ret)
}