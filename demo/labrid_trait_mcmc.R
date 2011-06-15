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

nchains <- 16
MaxTime = 1e6 
spec = list(alpha="global", sigma="indep", theta="global")

sfInit(parallel=T, cpu=16)
sfLibrary(wrightscape)
sfExportAll()
o <- sfLapply(1:nchains, function(i){ 
o <- phylo_mcmc(labrid$data['prot.y'], labrid$tree, pharyngeal,
                MaxTime=MaxTime, model_spec=spec, stepsizes=0.05)[[1]]
    })

burnin <- 1:1e3
chains <- o[[1]][-burnin,]
for(i in 2:nchains)
  chains <- rbind(chains, o[[i]][-burnin, ])

png(file="parameter_mcmc.png", width=3*480)
plot.phylo_mcmc(chains, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
dev.off()

upload("parameter_mcmc.png", script, gitaddr=gitaddr, tags=tags)

