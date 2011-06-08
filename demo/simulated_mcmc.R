require(wrightscape)
require(socialR)
source("parrotfish_data.R")
# since we can't install package while other reps are running
source("../R/wrightscape.R")
source("../R/mcmc.R")
source("../R/likelihood.R")
source("../R/prior_library.R")
# Create some simulated data on the parrotfish tree
true <- multiTypeOU(data=labrid$data["close"], tree=labrid$tree,regimes=intramandibular, 
                    model_spec=list(alpha="indep", sigma="global", theta="global"), 
                    Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                    method ="SANN", control=list(maxit=80000,temp=50,tmax=50))
true$Xo <- 1
true$alpha <- c(.01, 20)
true$sigma <- 10
true$theta <- 0
sim_trait <- simulate(true, seed=1)
## Start with a simple fit of indep alphas model to get some parameters
fit <- multiTypeOU(data=sim_trait, tree=labrid$tree, regimes=intramandibular, 
                  model_spec=list(alpha="indep", sigma="global", theta="global"), 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=80000,temp=50,tmax=50)) 



sfInit(parallel=T, cpu=4)
sfLibrary(wrightscape)
sfExportAll()
# MCMCMC the rc model
o <- phylo_mcmc(sim_trait, labrid$tree, intramandibular, MaxTime=1e5, indep=1e2,
                model_spec=list(alpha="indep", sigma="global", theta="global"),
                nchains=4)

png("alphadist.png")
cold_chain <- o$chains[[1]]
colnames(cold_chain) <- c("Pi", "Xo", "alpha1", "alpha2", "sigma", "theta")
poste_alpha1 <- density(cold_chain[, "alpha1"])
poste_alpha2 <- density(cold_chain[, "alpha2"])
xlim <- c(min(poste_alpha1$x, poste_alpha2$x), max(poste_alpha1$x, poste_alpha2$x)) 
ylim <- c(min(poste_alpha1$y, poste_alpha2$y), max(poste_alpha1$y, poste_alpha2$y)) 
plot(poste_alpha2, xlab="alpha", main="Selection Strength", xlim=xlim, ylim=ylim)
polygon(poste_alpha1, col=rgb(0,1,0,.5))
polygon(poste_alpha2, col=rgb(0,0,1,.5))
dev.off()

social_report(file="alphadist.png", tag="phylogenetics")


