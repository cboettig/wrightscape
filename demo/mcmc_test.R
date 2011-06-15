# mcmc_demo.R
rm(list=ls())
require(wrightscape)

##############
require(socialR)
script <- "mcmc_test.R"

tags <- c("phylogenetics parrotfish")
gitopts <- list(user = "cboettig", dir = "demo", repo = "wrightscape") 
gitaddr <- gitcommit(script, gitopts)
on.exit(system("git push")) #  For git links.  May prompt for pw,
tweet_errors(script, gitopts, tags)  ## tweet on error
#################

source("../R/mcmc.R")
source("../R/likelihood.R")
source("parrotfish_data.R")

MaxTime = 1e3 # 1e7 too great to store in mem, better start writing to file!
spec = list(alpha="global", sigma="indep", theta="global")

## uses [[1]] to return chains only, doesn't return the myCall
o <- phylo_mcmc(labrid$data['prot.y'], labrid$tree, intramandibular,
                MaxTime=MaxTime, model_spec=spec, stepsizes=0.05)[[1]]

plot.phylo_mcmc <- function(par_dist, ...){

 posterior <- vector("list", dim(par_dist)[2])
 subplot <- function(parname, ...){
  id <- grep(parname, colnames(par_dist))
  for(i in id){
    posterior[[i]] <- density(par_dist[,i])
  }
  xlim <- c(min(sapply(posterior[id], function(P) P$x)),
            max(sapply(posterior[id], function(P) P$x)))
  plot(posterior[[id[1]]], xlab=parname, xlim=xlim,
       )
  for(i in id)
    polygon(posterior[[id[1]]], col=c(0,1,0,.5))
 }

 par(mfrow=c(1,3))
 subplot("alpha", main="Selection Strength", ...)
 subplot("sigma", main="Trait Diversification", ...)
 subplot("theta", main="Optimum", ...)
}


png(file="parameter_mcmc.png", width=3*480)
plot.phylo_mcmc(o, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
dev.off()

upload("parameter_mcmc.png", script, gitaddr=gitaddr, tags=tags)


