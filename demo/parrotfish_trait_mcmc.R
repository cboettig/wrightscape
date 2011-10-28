# mcmc_demo.R
rm(list=ls())
require(wrightscape)
require(snowfall)

source("parrotfish_data.R")
MaxTime = 1e6
spec = list(alpha="fixed", sigma="indep", theta="global")
trait <-  "open"
nchains <- 1
burnin <- 1:1e3


# START SMART PLEASE
start <- multiTypeOU(data=labrid$data[trait], tree=labrid$tree, 
regimes=intramandibular, model_spec=spec) #,
#                  method ="SANN", control=list(maxit=100000,temp=50,tmax=20))


#sfInit(par=T, cpu=nchains)
sfLibrary(wrightscape)
sfExportAll()


o <- sfLapply(1:nchains, function(i){
chains <- phylo_mcmc(labrid$data[trait], labrid$tree, intramandibular,
		 MaxTime=MaxTime, model_spec=spec, stepsizes=0.05,
		 Xo=start$Xo, alpha=start$alpha, sigma=start$sigma,
		 theta=start$theta)
# returns [[1]]: chains, [[2]]: myCall, [[3]] colnames (for txtfile version)
chains[[1]]
})

# Concatinate chains 
chains <- o[[1]][-burnin,] # the first chain
for(i in 2:nchains)
	chains <- rbind(chains, o[[i]][-burnin, ])

# chains <- chains[[1]][-burnin,]  ## without parrallel

save(list=ls(), file="parrotfish_mcmc.Rdat")

  png(file="parameter_mcmc.png", width=3*480)
  plot.phylo_mcmc(chains, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
  dev.off()

#  upload("parameter_mcmc.png", script, gitaddr=gitaddr, 
#          tags=tags, comment=trait)

