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

boot_wright <- function(nboot=20){
	nboot = 1000
	require(phyloniche)
	require(wrightscape)
	data(bimac)
	tree <- with(bimac,ouchtree(nodes=node,ancestors=ancestor,times=time/max(time),labels=species))
	bm <- brown(log(bimac['size']), tree)
	ou1 <- hansen(log(bimac['size']), tree, bimac['OU.1'], 1, 1)


	lik_w <- model_bootstrap(ou1, NULL, nboot, update="wrightscape")


	model_list <- list(bm = bm, ou1 = ou1)
	LR_w <- LR_bootstrap_all(model_list, NULL, nboot=nboot, cpu=1, update="wrightscape")
	p_val_w <- t(summary(LR_w))

	save(list=ls(), file="corrected_ou.Rdat")
}

wright_plot <- function(LR_w){
	p_val_w <- t(summary(LR_w))
	pdf("BMvOU1.pdf")
	hist(-2*LR_w[[2]]$t, col="lightblue", border="white", xlab="Likelihood Ratio", cex.axis=1.6, cex.lab=1.6, main="")
	abline(v=-2*LR_w[[2]]$t0, lwd=4, lty=2, col="darkblue")
	text(-2*LR_w[[2]]$t0, 0.5*par()$yaxp[2], paste("p = ", round(p_val_w[2], digits=3)), cex=1.6)   # halfway up the vert line
	dev.off()

	pdf("OU1vBM.pdf")
	hist(-2*LR_w[[3]]$t, col="lightblue", border="white", xlab="Likelihood Ratio", cex.axis=1.6, cex.lab=1.6, main="")
	abline(v=-2*LR_w[[3]]$t0, lwd=4, lty=2, col="darkblue")
	text(-2*LR_w[[3]]$t0, 0.5*par()$yaxp[2], paste("p = ", round(p_val_w[3], digits=3)), cex=1.6)   # halfway up the vert line
	dev.off()
}
