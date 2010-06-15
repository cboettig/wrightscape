#traitmodels.R



traitmodels <- function(traits, tree){
	pars <- c(.44, .19, 2.9)
	fitpars <- c(1,1,0)


	traits[is.na(traits)] = 0 

	traits[1] = pars[3];

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
		as.double(traits),
		as.integer(states),
		as.integer(nstates),
		as.integer(n),
		as.double(pars),
		as.integer(fitpars),
		as.integer(npars),
		double(1)
	  )
	o[[10]]
}

testme <- function(){
	require(geiger)
	data(bimac)
	traits <- log(bimac$size)
	tree <- with(bimac,ouchtree(nodes=node,ancestors=ancestor,times=time,labels=species))
	traitmodels(traits, tree)
}
