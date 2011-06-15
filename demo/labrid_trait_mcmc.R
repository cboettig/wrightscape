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

comment=""

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

colnames(chains) <- c("Pi", "Xo", "alpha1", "alpha2", "sigma1", "sigma2", "theta")
par_dist <- chains



png(file="parameter_mcmc.png", width=3*480)
par(mfrow=c(1,3))
poste_alpha1 <- density(par_dist[, "alpha1"])
poste_alpha2 <- density(par_dist[, "alpha2"])
xlim <- c(min(poste_alpha1$x, poste_alpha2$x), max(poste_alpha1$x, poste_alpha2$x)) 
ylim <- c(min(poste_alpha1$y, poste_alpha2$y), max(poste_alpha1$y, poste_alpha2$y)) 
plot(poste_alpha2, xlab="alpha", main="Selection Strength", xlim=xlim, ylim=ylim, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
polygon(poste_alpha1, col=rgb(0,1,0,.5))
polygon(poste_alpha2, col=rgb(0,0,1,.5))

poste_theta1 <- density(par_dist[, "theta"])
plot(poste_theta1, xlab="theta", main="Optimum", cex=3, cex.lab=3, cex.main=3, cex.axis=3)
polygon(poste_theta1, col=rgb(0,1,0,.5))

poste_sigma1 <- density(par_dist[, "sigma1"])
poste_sigma1 <- density(par_dist[, "sigma2"])
xlim <- c(min(poste_sigma1$x, poste_sigma2$x), max(poste_sigma1$x, poste_sigma2$x)) 
plot(poste_sigma1, xlab="sigma", main="Diversification rate", cex=3, cex.lab=3, cex.main=3, cex.axis=3, xlim=xlim)
polygon(poste_sigma1, col=rgb(0,1,0,.5))
plot(poste_sigma2, xlab="sigma", main="Diversification rate", cex=3, cex.lab=3, cex.main=3, cex.axis=3, xlim=xlim)
polygon(poste_sigma2, col=rgb(0,0,1,.5))
dev.off()


upload("parameter_mcmc.png", script, gitaddr=gitaddr, tags=tags, comment=comment)

