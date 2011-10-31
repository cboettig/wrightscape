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
brownie <- multiTypeOU(data=labrid$data[trait], tree=labrid$tree, 
regimes=intramandibular, model_spec=spec) #,
#                  method ="SANN", control=list(maxit=100000,temp=50,tmax=20))


#sfInit(par=T, cpu=nchains)
sfLibrary(wrightscape)
sfExportAll()
out <- sfLapply(1:80, function(i){
	dat <- simulate(brownie) 
	out <- update(brownie, dat)
	out$alpha
})

