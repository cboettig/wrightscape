require(wrightscape)
source("parrotfish_data.R")
general = list(alpha="indep", sigma="indep", theta="indep")
rc = list(alpha="indep", sigma="global", theta="global")
fit_input <- list(data=labrid$data["close"], tree=labrid$tree,
                  regimes=intramandibular, model_spec=rc, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50))
true <- do.call(multiTypeOU, fit_input)
true$alpha <- c(.01, 20)
true$sigma <- 10
true$theta <- 0
sim_trait <- simulate(true)
## Now we're up and running with fake data

fit_input <- list(data=sim_trait, tree=labrid$tree,
                  regimes=intramandibular, model_spec=rc, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50))
fit <- do.call(multiTypeOU, sim_trait)
