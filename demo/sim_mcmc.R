require(wrightscape)
require(socialR)

script <- "sim_mcmc.R"
gitaddr <- gitcommit(script)
gitopts <- list(user = "cboettig", dir = "demo", repo = "wrightscape") 
on.exit(system("git push")) #  For git links.  May prompt for pw,
tags <- "phylogenetics"  ## multiple possible: space, delim, multiple items, etc.  
tweet_errors(script, gitopts, tags)  ## tweet on error


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
                model_spec=list(alpha="indep", sigma="global", theta="global"), 
                Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50)) 

nchains<-16
sfInit(parallel=T, cpu=8)
sfLibrary(wrightscape)
sfExportAll()
# MCMCMC the rc model

o <- sfLapply(1:nchains, function(i){ 
    phylo_mcmc(sim_trait, labrid$tree, intramandibular, 
               MaxTime=1e5, alpha=fit$alpha, sigma=fit$sigma, 
               theta=fit$theta, Xo=fit$Xo,
               model_spec=list(alpha="indep", sigma="global", theta="global"),
               stepsizes=0.5)[[1]]
    })

burnin <- 1:1e3
chains <- o[[1]][-burnin,]
for(i in 2:nchains)
  chains <- rbind(chains, o[[i]][-burnin, ])


png(file="parameter_mcmc.png", width=3*480)
plot.phylo_mcmc(chains, cex=3, cex.lab=3, cex.main=3, cex.axis=3)
dev.off()
upload("parameter_mcmc.png", script, gitaddr=gitaddr, tags=tags)


