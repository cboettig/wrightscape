# mcmc_demo.R
rm(list=ls())
require(wrightscape)

source("../R/mcmc.R")
source("../R/likelihood.R")
source("parrotfish_data.R")

MaxTime = 1e3 # 1e7 too great to store in mem, better start writing to file!
spec = list(alpha="global", sigma="indep", theta="global")

## uses [[1]] to return chains only, doesn't return the myCall
o <- phylo_mcmc(labrid$data['prot.y'], labrid$tree, intramandibular,
                MaxTime=MaxTime, model_spec=spec, stepsizes=0.05)[[1]]

plot.phylo_mcmc(o, cex=3, cex.lab=3, cex.main=3, cex.axis=3)



