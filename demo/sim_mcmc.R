rm(list=ls())
require(wrightscape)
require(socialR)

script <- "sim_mcmc.R"
gitaddr <- gitcommit(script)
gitopts <- list(user = "cboettig", dir = "demo", repo = "wrightscape") 
on.exit(system("git push")) #  For git links.  May prompt for pw,
tags <- "phylogenetics"  ## multiple possible: space, delim, multiple items, etc.  
tweet_errors(script, gitopts, tags)  ## tweet on error


source("parrotfish_data.R")
# since we can't install package while other reps are running

# Create some simulated data on the parrotfish tree
true <- multiTypeOU(
        data=labrid$data["close"], tree=labrid$tree,regimes=intramandibular, 
        model_spec=list(alpha="indep", sigma="global", theta="global"), 
        Xo=NULL, alpha = .1, sigma = .1, theta=NULL)
true$Xo <- 1
true$alpha <- c(.01, 20)
true$sigma <- 10
true$theta <- 0
sim_trait <- simulate(true, seed=1)

#nchains<-8
#sfInit(parallel=T, cpu=8)
#sfLibrary(wrightscape)
#sfExportAll()
# MCMCMC the rc model

#o <- sfLapply(1:nchains, function(i){ 
 chains<-   phylo_mcmc(sim_trait, labrid$tree, intramandibular, 
               MaxTime=1e5, alpha=true$alpha, sigma=true$sigma, 
               theta=true$theta, Xo=true$Xo,
               model_spec=list(alpha="indep", sigma="global", theta="global"),
               stepsizes=0.05)[[1]]
   # })

#burnin <- 1:1e3
#chains <- o[[1]][-burnin,]
#for(i in 2:nchains)
#  chains <- rbind(chains, o[[i]][-burnin, ])


png(file="parameter_mcmc.png", width=3*480)
plot.phylo_mcmc(chains, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
dev.off()
upload("parameter_mcmc.png", script, gitaddr=gitaddr, tags=tags)


