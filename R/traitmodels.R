#traitmodels.R
traitmodels <- function(){

	data(bimac)
	traits <- log(bimac$size)
	tree <- with(bimac,ouchtree(nodes=node,ancestors=ancestor,times=time,labels=species))


#	if(class(tree) == "phylo"){
#		tree = ape2ouch(tree);
#	}
#	if(class(data) == "data.frame"){
#		traits = data[,1]
#	} else {
#		traits = as.numeric(c( rep('NA',length(data)-1)  , data));
#	}

	traits[is.na(traits)] = 0 

	## Make this into a function that takes data from an OU tree, optionally with a painting, and a mapping of that painting to models
#	if(class(tree) == "ouch"){
		ancestors <- as.numeric(tree@ancestors)
		ancestors[is.na(ancestors)] = 0 

		plot(tree)

		times <-tree@times
		## ouch gives cumulative time, not branch-length!!
		t <- times
		for(i in 1:length(ancestors) ){
			if(ancestors[i] > 0){
				times[i] <- times[i] - t[ancestors[i]]
			}
		}
		ancestors <- ancestors-1
		n = length(times)
	
#	}



	.C("traitmodels",
		as.numeric(times),
		as.integer(ancestors),
		as.numeric(traits),
		as.integer(states),
		as.integer(n),
		as.numeric(pars),
		as.integer(fitpars),
		as.integer(npars)
	  )
}
