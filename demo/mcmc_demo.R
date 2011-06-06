# mcmc_demo.R
require(wrightscape)
require(socialR)
tags <- c("phylogenetics mcmcmc parrotfish")
source("parrotfish_data.R")
sfInit(parallel=T, cpu=4)
sfLibrary(wrightscape)
sfExportAll()

o <- general_mcmc(labrid$data['gape.y'], labrid$tree, intramandibular,
                  alpha=.1, sigma=.1, MaxTime=1e6, indep=1e2)


cold_chain <- o[[1]]
colnames(cold_chain) <- c("Pi", "Xo", "alpha1", "alpha2", "sigma1", "sigma2", "theta1", "theta2")


png("posteriors.png", width=3*480)
par(mfrow=c(1,3))
poste_alpha1 <- density(cold_chain[, "alpha1"])
poste_alpha2 <- density(cold_chain[, "alpha2"])
xlim <- c(min(poste_alpha1$x, poste_alpha2$x), max(poste_alpha1$x, poste_alpha2$x)) 
ylim <- c(min(poste_alpha1$y, poste_alpha2$y), max(poste_alpha1$y, poste_alpha2$y)) 
plot(poste_alpha2, xlab="alpha", main="Selection Strength", xlim=xlim, ylim=ylim)
polygon(poste_alpha1, col=rgb(0,1,0,.5))
polygon(poste_alpha2, col=rgb(0,0,1,.5))

poste_theta1 <- density(cold_chain[, "theta1"])
poste_theta2 <- density(cold_chain[, "theta2"])
xlim <- c(min(poste_theta1$x, poste_theta2$x), max(poste_theta1$x, poste_theta2$x))
ylim <- c(min(poste_theta1$y, poste_theta2$y), max(poste_theta1$y, poste_theta2$y)) 
plot(poste_theta2, xlab="theta", main="Optimum", xlim=xlim, ylim=ylim)
polygon(poste_theta1, col=rgb(0,1,0,.5))
polygon(poste_theta2, col=rgb(0,0,1,.5))

poste_sigma1 <- density(cold_chain[, "sigma1"])
poste_sigma2 <- density(cold_chain[, "sigma2"])
xlim <- c(min(poste_sigma1$x, poste_sigma2$x), max(poste_sigma1$x, poste_sigma2$x)) 
ylim <- c(min(poste_sigma1$y, poste_sigma2$y), max(poste_sigma1$y, poste_sigma2$y)) 
plot(poste_sigma2, xlab="sigma", main="Diversification rate", xlim=xlim, ylim=ylim)
polygon(poste_sigma1, col=rgb(0,1,0,.5))
polygon(poste_sigma2, col=rgb(0,0,1,.5))
dev.off()

social_report(file="posteriors.png", tag=tags)


png("convergenceTemp.png")
  burnin <- 1:1e5
  plot(o[[4]][-burnin,1], type="l", col=rgb(1,0,0.3))
  lines(o[[3]][-burnin,1], col=rgb(1,0,0.5))
  lines(o[[2]][-burnin,1], col=rgb(1,0,0.8))
  lines(o[[1]][-burnin,1])
dev.off()

social_report(file="convergenceTemp.png", tag=tags, comment="MaxTime=1e5, indep=1e1, stepsizes=.2")



