# File: mcmc.R
# Author: Carl Boettiger <cboettig@gmail.com>
# Date: 2011-11-16
# License: BSD
#
# These functions are stull under development, not yet exported to user

phylo_mcmc <- function(data, tree, regimes, model_spec =
                       list(alpha="indep", sigma="indep", theta="indep"),
                       Xo=NULL, alpha=1, sigma=1, theta=NULL, prior=NULL,
                       MaxTime, stepsizes=.1, write=NULL){

  myCall <- match.call() # keep input for the record
  
  n_regimes <- length(levels(regimes))
  if(is.null(prior))
    prior <- flat_prior_creator(model_spec, n_regimes)
  par <- setup_pars(data, tree, regimes, model_spec, Xo=Xo, 
                    alpha=alpha, sigma=sigma, theta=theta)
  f <- llik.closure(data, tree, regimes, model_spec)
  chain <- mcmc_fn(par, f, prior, MaxTime=MaxTime, 
                   stepsizes=stepsizes, write=write)


  heading <- unique_names(model_spec, n_regimes) 
  if(is.null(write))
    colnames(chain) <- heading 
  else
    chain=write
  list(chain=chain, myCall=myCall, colnames=heading)
}

unique_names <- function(model_spec, n_regimes){
  indices <- get_indices(model_spec, n_regimes)
  indices <- lapply(indices, function(x) x[!duplicated(x)] )
  tmp <- c("Pi" = 1, "Xo"=2, "alpha" = indices$alpha_i,
           "sigma"=indices$sigma_i, "theta"=indices$theta_i)
  names(tmp)
}

phylo_mcmcmc <- function(data, tree, regimes, model_spec =
                       list(alpha="indep", sigma="indep", theta="indep"),
                       Xo=NULL, alpha=1, sigma=1, 
                       theta=NULL, prior=NULL, MaxTime, 
                       indep=100, stepsizes=.1, nchains=4, 
                       Delta_T =1, ...){

  myCall <- match.call() # keep input for the record
  
  n_regimes <- length(levels(regimes))
  if(is.null(prior))
    prior <- flat_prior_creator(model_spec, n_regimes)
  par <- setup_pars(data, tree, regimes, model_spec, Xo=Xo, 
                    alpha=alpha, sigma=sigma, theta=theta)
  f <- llik.closure(data, tree, regimes, model_spec)

  ## Assemble starting points of the different chains
  ## Could be much more creative/clever than this
  if(nchains > 1){
  pars <- c(par, lapply(2:nchains, 
                 function(i){
                    sapply(rexp(length(par), 5), abs)
                 }))
  } else {
    pars <- list(par)
  }
  chains <- mcmcmc_fn(pars, f, prior, MaxTime=MaxTime, indep=indep,
                      stepsizes=stepsizes, Delta_T = Delta_T, ...) 
  list(chains=chains, myCall=myCall)
}






## Verify that alpha1, sigma1, etc corresponds to regimes[1], etc!!!
## make sure that colnames are given in increasing numerical order!
plot.phylo_mcmc <- function(par_dist, xlim=NULL, ...){
 posterior <- vector("list", dim(par_dist)[2])

 colors <- c( rgb(0,1,0.5), rgb(0,0,1,.5), rgb(1,0,0,.5), rgb(1,1,0.5) )
 
 subplot <- function(parname, ...){
  id <- grep(parname, colnames(par_dist))
  for(i in id){
    posterior[[i]] <- density(par_dist[,i])
  }
  if(is.null(xlim))
    xlim <- c(min(sapply(posterior[id], function(P) P$x)),
              max(sapply(posterior[id], function(P) P$x)))
  plot(posterior[[id[1]]], xlab=parname, xlim=xlim, ...)
  cindex <- 1
  for(i in id){
    polygon(posterior[[i]], col=colors[cindex])
    cindex <- cindex + 1
  }
 }

 par(mfrow=c(1,3))
 subplot("alpha", main="Selection Strength", ...)
 subplot("sigma", main="Trait Diversification", ...)
 subplot("theta", main="Optimum", ...)
}





