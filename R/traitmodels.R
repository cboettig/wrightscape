#traitmodels.R



traitmodels <- function(tree, data){

	pars <- c(tree@sigma, tree@sqrt.alpha, tree@theta[[1]][1] )
	fitpars <- c(1,1,1)

	data <- data[[1]]
	data[is.na(data)] = 0 


	ancestors <- as.numeric(tree@ancestors)
	ancestors[is.na(ancestors)] = 0 

	## ouch gives cumulative time, not branch-length!!
	anc <- as.integer(tree@ancestors[!is.na(tree@ancestors)])
	lengths <- c(0, tree@times[!is.na(tree@ancestors)] - tree@times[anc] )
	times <- lengths/max(tree@times)

	ancestors <- ancestors-1
	n <- length(times)

	states <- rep(0, n)
	nstates <- 1
	npars <- length(pars)
		
	o<- .C("traitmodels",
		as.double(times),
		as.integer(ancestors),
		as.double(data),
		as.integer(states),
		as.integer(nstates),
		as.integer(n),
		as.double(pars),
		as.integer(fitpars),
		as.integer(npars),
		double(1)
	  )
	list(loglik=o[[10]], sigma = sqrt( o[[7]][1] ), sqrt.alpha = sqrt( o[[7]][2] ), theta = sqrt( o[[7]][3] ))  
}

testme <- function(){
	require(phyloniche)
	data(bimac)
	tree <- with(bimac,ouchtree(nodes=node,ancestors=ancestor,times=time/max(time),labels=species))
	bm <- brown(log(bimac['size']), tree)
	ou1 <- hansen(log(bimac['size']), tree, bimac['OU.1'], 1, 1)
	model_list <- list(bm = bm, ou1 = ou1)
	o <- LR_bootstrap_all(model_list, NULL, nboot=20, cpu=2, update="wrightscape")
}



