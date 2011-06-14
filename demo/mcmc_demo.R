# mcmc_demo.R
require(wrightscape)
require(socialR)
tags <- c("phylogenetics parrotfish")
script <- "mcmc_demo.R"
gitcommit(script)
gitopts = list(user = "cboettig", dir = "demo", repo = "wrightscape") 
on.exit(system("git push")) #  For git links.  May prompt for pw,
tweet_errors(script, gitopts, tags)  ## tweet on error

source("parrotfish_data.R")
MaxTime = 1e6 # 1e7 too great to store in mem, better start writing to file!
spec = list(alpha="indep", sigma="global", theta="global")

comment="open"

sfInit(parallel=T, cpu=8)
sfLibrary(wrightscape)
sfExportAll()
o <- sfLapply(1:nchains, function(i){ 
o <- phylo_mcmc(labrid$data['prot.y'], labrid$tree, intramandibular,
                MaxTime=MaxTime, model_spec=spec, stepsizes=0.05)[[1]]
    })

burnin <- 1:1e3
chains <- o[[1]][-burnin,]
for(i in 2:nchains)
  chains <- rbind(chains, o[[i]][-burnin, ])

colnames(chains) <- c("Pi", "Xo", "alpha1", "alpha2", "sigma", "theta")
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

poste_sigma1 <- density(par_dist[, "sigma"])
plot(poste_sigma1, xlab="sigma", main="Diversification rate", cex=3, cex.lab=3, cex.main=3, cex.axis=3)
polygon(poste_sigma1, col=rgb(0,1,0,.5))
dev.off()


upload("parameter_mcmc.png", script, gitopts=gitopts, tags=tags, comment=comment)

