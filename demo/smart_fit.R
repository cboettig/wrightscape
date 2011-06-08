require(wrightscape)
require(socialR)
source("parrotfish_data.R")

# since we can't install package while other reps are running
source("../R/wrightscape.R")
source("../R/mcmc.R")
source("../R/likelihood.R")
source("../R/prior_library.R")

brownie = list(alpha="fixed", sigma="indep", theta="globa")
rc = list(alpha="indep", sigma="global", theta="global")
fit_input <- list(data=labrid$data["close"], tree=labrid$tree,
                  regimes=intramandibular, model_spec=rc, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50))
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
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL) 

gen_fit <- smart_multiType(data=sim_trait, tree=labrid$tree,
                  regimes=intramandibular, model_spec=brownie, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL)
require(pmc)
sfInit(parallel=T, cpu=16)
sfLibrary(wrightscape)
sfExportAll()
boots <- montecarlotest(fit, gen_fit, nboot=80, cpu=16, GetParNames=TRUE)
social_plot(plot(boots), tags="phylogenetics")


