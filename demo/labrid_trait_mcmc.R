# mcmc_demo.R
rm(list=ls())
require(wrightscape)

##############
require(socialR)
script <- "labrid_trait_mcmc.R"
tags <- c("phylogenetics parrotfish")
gitopts <- list(user = "cboettig", dir = "demo", repo = "wrightscape") 
gitaddr <- gitcommit(script, gitopts)
on.exit(system("git push")) #  For git links.  May prompt for pw,
tweet_errors(script, gitopts, tags)  ## tweet on error
#################

source("labrid_data.R")

MaxTime = 1e6 
spec = list(alpha="indep", sigma="indep", theta="global")
traits <- c("prot.y", "close", "open", "gape.y")

nchains <- 4
burnin <- 1:1e5


sfInit(parallel=T, cpu=4)
sfLibrary(wrightscape)
sfLibrary(socialR)
sfExportAll()

#sfLapply(traits, function(trait){
  # START SMART PLEASE
  start <- multiTypeOU(data=labrid$data[trait], tree=labrid$tree, 
                  regimes=pharyngeal, model_spec=spec,
                  method ="SANN", control=list(maxit=100000,temp=50,tmax=20))

  o <- sfLapply(1:4, function(i)){
    chains <- phylo_mcmc(labrid$data[trait], labrid$tree, pharyngeal,
                         MaxTime=MaxTime, model_spec=spec, stepsizes=0.05,
                         Xo=start$Xo, alpha=start$alpha, sigma=start$sigma,
                         theta=start$theta)
    # returns [[1]]: chains, [[2]]: myCall, [[3]] colnames (for txtfile version)
  chains[[1]]
  }

  # the first chain
  chains <- o[[1]][-burnin,]
  for(i in 2:nchains){
    chains <- rbind(chains, o[[i]][-burnin, ]
  }

  # chains <- chains[[1]][-burnin,]  ## without parrallel


  png(file="parameter_mcmc.png", width=3*480)
  plot.phylo_mcmc(chains, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
  dev.off()
  upload("parameter_mcmc.png", script, gitaddr=gitaddr, 
          tags=tags, comment=trait)

})
