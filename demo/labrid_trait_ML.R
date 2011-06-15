# parrotfish.R
require(wrightscape)
require(pmc)

############ Notebook logging header ##############
require(socialR)
script <- "labrid_trait_ML.R"
tags="phylogenetics"
gitopts <- list(user = "cboettig", dir = "demo", repo = "wrightscape") 
gitaddr <- gitcommit(script, gitopts)
on.exit(system("git push")) 
tweet_errors(script, gitopts, tags) 
####################################################

source("labrid_data.R")

alphas <- multiTypeOU(data=labrid$data["prot.y"], tree=labrid$tree, regimes=pharyngeal, 
            model_spec=list(alpha="indep", sigma="global", theta="indep"), 
            Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
            method ="SANN", control=list(maxit=50000,temp=50,tmax=20))

sigmas <- multiTypeOU(data=labrid$data["prot.y"], tree=labrid$tree, regimes=pharngeal, 
                model_spec=list(alpha="global", sigma="indep", theta="indep"), 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=50000,temp=50,tmax=20))



require(snowfall)
sfInit(parallel=TRUE, cpu=8)
sfLibrary(wrightscape)
sfExportAll()

boots <- montecarlotest(sigmas, alphas, nboot=400, cpu=8)
png("sigmas_v_alphas.png")
plot(boots)
dev.off()

upload("sigmas_v_alphas.png", script=script, tags=tags, gitaddr=gitaddr)

## this needs to become a smart function
finalplots <- function(boots){
  par_dist <- t(boots$test_par_dist) 

 png("bootstrap_pars.png", width=3*480) 
    par(mfrow=c(1,3))
    poste_alpha1 <- density(par_dist[, "alpha1"])
    poste_alpha2 <- density(par_dist[, "alpha2"])
    xlim <- c(min(poste_alpha1$x, poste_alpha2$x), max(poste_alpha1$x, poste_alpha2$x)) 
    ylim <- c(min(poste_alpha1$y, poste_alpha2$y), max(poste_alpha1$y, poste_alpha2$y)) 
    plot(poste_alpha2, xlab="alpha", main="Selection Strength", xlim=xlim, ylim=ylim, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
    polygon(poste_alpha1, col=rgb(0,1,0,.5))
    polygon(poste_alpha2, col=rgb(0,0,1,.5))

    poste_theta1 <- density(par_dist[, "theta1"])
    poste_theta2 <- density(par_dist[, "theta2"])
    xlim <- c(min(poste_theta1$x, poste_theta2$x), max(poste_theta1$x, poste_theta2$x))
    ylim <- c(min(poste_theta1$y, poste_theta2$y), max(poste_theta1$y, poste_theta2$y)) 
    plot(poste_theta2, xlab="theta", main="Optimum", xlim=xlim, ylim=ylim, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
    polygon(poste_theta1, col=rgb(0,1,0,.5))
    polygon(poste_theta2, col=rgb(0,0,1,.5))

    poste_sigma1 <- density(par_dist[, "sigma1"])
    poste_sigma2 <- density(par_dist[, "sigma2"])
    xlim <- c(min(poste_sigma1$x, poste_sigma2$x), max(poste_sigma1$x, poste_sigma2$x)) 
    ylim <- c(min(poste_sigma1$y, poste_sigma2$y), max(poste_sigma1$y, poste_sigma2$y)) 
    plot(poste_sigma2, xlab="sigma", main="Diversification rate", xlim=xlim, ylim=ylim, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
    polygon(poste_sigma1, col=rgb(0,1,0,.5))
    polygon(poste_sigma2, col=rgb(0,0,1,.5))
  dev.off()

upload("bootstrap_pars.png", script=script, gitaddr=gitaddr, tags=tags)

  par_dist <- t(boots$null_par_dist) 
  png(file="parameter_bootstraps.png", width=3*480)
    par(mfrow=c(1,3))
    poste_alpha1 <- density(par_dist[, "alpha1"])
    poste_alpha2 <- density(par_dist[, "alpha2"])
    xlim <- c(min(poste_alpha1$x, poste_alpha2$x), max(poste_alpha1$x, poste_alpha2$x)) 
    ylim <- c(min(poste_alpha1$y, poste_alpha2$y), max(poste_alpha1$y, poste_alpha2$y)) 
    plot(poste_alpha2, xlab="alpha", main="Selection Strength", xlim=xlim, ylim=ylim, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
    polygon(poste_alpha1, col=rgb(0,1,0,.5))
    polygon(poste_alpha2, col=rgb(0,0,1,.5))

    poste_theta1 <- density(par_dist[, "theta1"])
    poste_theta2 <- density(par_dist[, "theta2"])
    xlim <- c(min(poste_theta1$x, poste_theta2$x), max(poste_theta1$x, poste_theta2$x))
    ylim <- c(min(poste_theta1$y, poste_theta2$y), max(poste_theta1$y, poste_theta2$y)) 
    plot(poste_theta2, xlab="theta", main="Optimum", xlim=xlim, ylim=ylim, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
    polygon(poste_theta1, col=rgb(0,1,0,.5))
    polygon(poste_theta2, col=rgb(0,0,1,.5))

    poste_sigma1 <- density(par_dist[, "sigma1"])
    poste_sigma2 <- density(par_dist[, "sigma2"])
    xlim <- c(min(poste_sigma1$x, poste_sigma2$x), max(poste_sigma1$x, poste_sigma2$x)) 
    ylim <- c(min(poste_sigma1$y, poste_sigma2$y), max(poste_sigma1$y, poste_sigma2$y)) 
    plot(poste_sigma2, xlab="sigma", main="Diversification rate", xlim=xlim, ylim=ylim, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
    polygon(poste_sigma1, col=rgb(0,1,0,.5))
    polygon(poste_sigma2, col=rgb(0,0,1,.5))
  dev.off()
  upload("parameter_bootstraps.png", script=script, tags=tags, gitaddr=gitaddr)
}

finalplots(boots)



