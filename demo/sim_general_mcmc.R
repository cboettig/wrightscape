require(wrightscape)
require(socialR)
source("parrotfish_data.R")
# since we can't install package while other reps are running
source("../R/wrightscape.R")
source("../R/mcmc.R")
source("../R/likelihood.R")
source("../R/prior_library.R")
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
## Start with a simple fit of indep alphas model to get some parameters
## of course we could also start with the true values.  
fit <- multiTypeOU(data=sim_trait, tree=labrid$tree, regimes=intramandibular, 
                model_spec=list(alpha="indep", sigma="indep", theta="indep"), 
                Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50)) 

nchains<-8
sfInit(parallel=T, cpu=8)
sfLibrary(wrightscape)
sfExportAll()
# MCMCMC the rc model

o <- sfLapply(1:nchains, function(i){ 
    phylo_mcmc(sim_trait, labrid$tree, intramandibular, 
               MaxTime=1e5, alpha=fit$alpha, sigma=fit$sigma, 
               theta=fit$theta, Xo=fit$Xo,
               model_spec=list(alpha="indep", sigma="indep", theta="indep"),
               stepsizes=0.5)[[1]]
    })

burnin <- 1:1e3
chains <- o[[1]][-burnin,]
for(i in 2:nchains)
  chains <- rbind(chains, o[[i]][-burnin, ])
colnames(chains) <- c("Pi", "Xo", "alpha1", "alpha2", "sigma1", "sigma2", "theta1", "theta2")
par_dist <- chains


run_id <- runif(1,1e9,1e10)
script_name <- paste("run_", run_id, ".txt", sep="")
savehistory(file=script_name)
script <- paste(readLines(script_name), collapse = "        ")




social_plot({
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
}, file="parameter_boostraps.png", width=3*480, tag="phylogenetics", comment=script)


