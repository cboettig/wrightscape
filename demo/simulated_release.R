require(wrightscape)
require(socialR)
source("parrotfish_data.R")

# since we can't install package while other reps are running
source("../R/wrightscape.R")
source("../R/mcmc.R")
source("../R/likelihood.R")
source("../R/prior_library.R")

brownie = list(alpha="fixed", sigma="indep", theta="global")
rc = list(alpha="indep", sigma="global", theta="global")
fit_input <- list(data=labrid$data["close"], tree=labrid$tree,
                  regimes=intramandibular, model_spec=rc, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL)
true <- do.call(multiTypeOU, fit_input)
true$Xo <- 1
true$alpha <- c(.01, 20)
true$sigma <- 10
true$theta <- 0
sim_trait <- simulate(true, seed=1)

## Repeat against a standard Nelder-Mead vs SANN
## Now we're up and running with fake data
fit <- multiTypeOU(data=sim_trait, tree=labrid$tree,
                  regimes=intramandibular, model_spec=rc, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL) #,
#                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50))

brownie_fit <- multiTypeOU(data=sim_trait, tree=labrid$tree,
                  regimes=intramandibular, model_spec=brownie, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL)
#,
 #                 method ="SANN", control=list(maxit=80000,temp=25,tmax=50))

print(getParameters.multiOU(brownie_fit))

#simulate(brownie_fit) -> X

#update(brownie_fit, simulate(brownie_fit))

#require(pmc)
#sfInit(parallel=F)
#sfLibrary(wrightscape)
#sfExportAll()
#boots <- montecarlotest(brownie_fit, fit, nboot=80, cpu=16, GetParNames=FALSE)
#social_plot(plot(boots), tags="phylogenetics")


####### Plot the distributions
function(boots){
rownames(boots$test_par_dist) <- names(getParameters(gen_fit))
par_dist <- t(boots$test_par_dist) 

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
}, file="parameter_boostraps.png", width=3*480, tag="phylogenetics")


rownames(boots$null_par_dist) <- names(getParameters(fit))
par_dist <- t(boots$null_par_dist)
social_plot({
poste_alpha1 <- density(par_dist[, "alpha1"])
poste_alpha2 <- density(par_dist[, "alpha2"])
xlim <- c(min(poste_alpha1$x, poste_alpha2$x), max(poste_alpha1$x, poste_alpha2$x)) 
ylim <- c(min(poste_alpha1$y, poste_alpha2$y), max(poste_alpha1$y, poste_alpha2$y)) 
plot(poste_alpha2, xlab="alpha", main="Selection Strength", xlim=xlim, ylim=ylim, cex=2, cex.lab=2, cex.main=2, cex.axis=2)
polygon(poste_alpha1, col=rgb(0,1,0,.5))
polygon(poste_alpha2, col=rgb(0,0,1,.5))
}, file="parameter_boostraps.png", width=480, tag="phylogenetics")




}










