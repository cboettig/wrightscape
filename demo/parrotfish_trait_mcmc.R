# mcmc_demo.R
rm(list=ls())
require(wrightscape)

##############
require(socialR)
script <- "parrotfish_trait_mcmc.R"
tags <- c("phylogenetics")
gitopts <- list(user = "cboettig", dir = "demo", repo = "wrightscape") 
gitaddr <- gitcommit(script, gitopts)
on.exit(system("git push")) #  For git links.  May prompt for pw,
tweet_errors(script, gitopts, tags)  ## tweet on error
#################

source("parrotfish_data.R")

MaxTime = 1e6 # 1e7 too great to store in mem, better start writing to file!
spec = list(alpha="indep", sigma="global", theta="global")
traits <- c("prot.y")# , "close", "open", "gape.y")



sfInit(parallel=T, cpu=4)
sfLibrary(wrightscape)
sfLibrary(socialR)
sfExportAll()

sfLapply(traits, function(trait){
  # START SMART PLEASE
  start <- multiTypeOU(data=labrid$data[trait], tree=labrid$tree, 
                  regimes=intramandibular, model_spec=spec,
                  method ="SANN", control=list(maxit=100000,temp=50,tmax=20))

  print(getParameters(start))

  chains <- phylo_mcmc(labrid$data[trait], labrid$tree, intramandibular,
                       MaxTime=MaxTime, model_spec=spec, stepsizes=0.05,
                       Xo=start$Xo, alpha=start$alpha, sigma=start$sigma,
                       theta=start$theta)
  # returns [[1]]: chains, [[2]]: myCall, [[3]] colnames (for txtfile version)

  burnin <- 1:1e5
  chains <- chains[[1]][-burnin,]
  png(file="parameter_mcmc.png", width=3*480)
  plot.phylo_mcmc(chains, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
  dev.off()
  upload("parameter_mcmc.png", script, gitaddr=gitaddr, 
          tags=tags, comment=trait)

})
