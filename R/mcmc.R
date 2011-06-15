######### CONSIDER IMPORTING THESE ALL FROM "mcmcTools" package instead ###########
step_fn <- function(pars, stepsizes = rep(.02, length(pars))){
# Sequential random updating 
  j <- sample(1:length(pars), 1)
  pars[j] <- rnorm(1, pars[j], stepsizes[j])
  pars
}

## Proposal density is symmetric, so we won't need Q
Q <- function(pars, proposed){
  dnorm(pars, proposed, stepsizes=1, log=TRUE)
#    alpha <- loglik(pars) + prior(pars) + Q(pars, proposed) - loglik(proposed) - prior(proposed) - Q(proposed, pars)
}

# The basic mcmc function
mcmc_fn <- function(pars, loglik, prior, MaxTime=1e3, stepsizes=.02, ...){
  if(length(stepsizes)==1)
    stepsizes <- rep(stepsizes, length(pars))
  history <- matrix(NA, nrow=MaxTime, ncol=(1+length(pars)))
  for(t in 1:MaxTime){
    Pi <- loglik(pars)+prior(pars)
    history[t,] <- c(Pi, pars)
    proposed <- step_fn(pars, stepsizes)
    alpha <- exp(loglik(proposed)+prior(proposed) - Pi)
  if (alpha > runif(1) )
      pars <- proposed
  }
  history
}

beta <- function(i, Delta_T=1){
  1/(1+Delta_T*(i-1))
}


mcmcmc_fn <- function(pars, loglik, prior, MaxTime=1e3, indep=100, stepsizes=.02, Delta_T=1, ...){
# Metropolis Coupled Markov Chain Monte Carlo
# Args:
#   pars: a list of length n_chains, with numerics pars[[i]] that can be passed to loglik
#   loglik: a function to calculate the log-likelihood of chain i at pars[[i]], 
#   prior: a function to calculate the prior density
#   MaxTime: length of time to run the chain
#   indep: period of time for which chains should wander independently
#   stepsizes: step of proposal distribution (can be numeric of length 1 or length pars)
#   Delta_T: amount heated chains are increased, 0 = all cold.  
# Returns:
#   chains: list containing matrix for each chain, first col is loglik + log prior prob,
#           remaining columns are fn parameters in order given in the pars[[i]]
  n_chains <- length(pars)
  n_pars <- length(pars[[1]])

  # in case we want to store the complete history.  Should have the option
  # of writing this to a file for speed? Or do that all in C...
  chains <- lapply(1:n_chains, function(i)  matrix(NA, nrow=MaxTime, ncol=(1+n_pars)) )

  # The independent intervals, lets us run chains in parallel during these periods
  Interval <- matrix(1:MaxTime, nrow=indep)

  # Note the outer time loop over intervals, and 
  # an inner time loop that can be parallelized over chains
  for(s in 1:(MaxTime/indep)){
  # Evolve chains independently for "indep" time steps
    out <- sfLapply(1:n_chains, 
            function(i){
              out <- matrix(NA, ncol=n_pars+1, nrow=indep)
              # Inner time loop
              for(t in 1:indep){
                Pi <- loglik(pars[[i]]) + prior(pars[[i]])
                out[t,] <- c(Pi, pars[[i]]) # more simply could print this to file, save mem
                proposed <- step_fn(pars[[i]], stepsizes)
                # Normal Hastings ratio weighted by temp fn beta 
                alpha <- exp( beta(i, Delta_T) * ( loglik(proposed)+prior(proposed) - Pi ) )
                if ( alpha  > runif(1) )
                  pars[[i]] <- proposed
              }
              out
            })
    # write to output
    for(i in 1:n_chains){
      chains[[i]][Interval[,s],] <- out[[i]]
      pars[[i]] <- out[[i]][indep,][-1] # copy parameters (but not loglik)
    }
    # This isn't quite precise bc this isn't being counted as a timestep yet!
    # Propose a chain swap every "indep" time steps beween pars_i and j
    pick <- sample(1:length(pars), 2)
    i <- pick[1]; j <- pick[2] # for convience
    R <- 
      beta(i, Delta_T) * ( loglik(pars[[j]]) + prior(pars[[j]]) - loglik(pars[[i]]) - prior(pars[[i]]) ) +
      beta(j, Delta_T) * ( loglik(pars[[i]]) + prior(pars[[i]])  - loglik(pars[[j]]) - prior(pars[[j]]) )
    ## verbose output about swaps
    # print(paste("swap chain", i, "with", j, "proposed, R =", exp(R)))
    if(exp(R) > runif(1)){
      # print("swap accepted")
      pars[[i]] <- pars[[j]] # accept the swap
    }

  }
  # Returns the full history of all chains
  chains 
}



########## The actual Phylogenetic Model MCMCMC #############
########## Depends on the above functions, could depend on mcmcTools instead ########

## Should take number of chains as an option
phylo_mcmc <- function(data, tree, regimes, model_spec =
                       list(alpha="indep", sigma="indep", theta="indep"),
                       Xo=NULL, alpha=1, sigma=1, 
                       theta=NULL, prior=NULL, MaxTime, 
                       indep=100, stepsizes=.1, 
                       Delta_T =1, ...){

  myCall <- match.call() # keep input for the record
  
  n_regimes <- length(levels(regimes))
  if(is.null(prior))
    prior <- flat_prior_creator(model_spec, n_regimes)
  par <- setup_pars(data, tree, regimes, model_spec, Xo=Xo, 
                    alpha=alpha, sigma=sigma, theta=theta)
  f <- llik.closure(data, tree, regimes, model_spec)
  chain <- mcmc_fn(par, f, prior, MaxTime=MaxTime, indep=indep,
                      stepsizes=stepsizes, Delta_T = Delta_T, ...)


  colnames(chain) <- unique_names(model_spec, n_regimes) 
  list(chain=chain, myCall=myCall)
}

unique_names <- function(model_spec, n_regimes){
  indices <- get_indices(model_spec, n_regimes)
  indices <- lapply(indices, function(x) x[!duplicated(x)] )
  tmp <- c("Pi" = 1, "Xo"=2, "alpha" = indices$alpha_i,
           "sigma"=indices$sigma_i, "theta"=indices$theta_i)
  names(tmp)
}

## Should take number of chains as an option
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
plot.phylo_mcmc <- function(par_dist, ...){
 posterior <- vector("list", dim(par_dist)[2])

 colors <- c( rgb(0,1,0.5), rgb(0,0,1,.5), rgb(1,0,0,.5), rgb(1,1,0.5) )
 
 subplot <- function(parname, ...){
  id <- grep(parname, colnames(par_dist))
  for(i in id){
    posterior[[i]] <- density(par_dist[,i])
  }
  xlim <- c(min(sapply(posterior[id], function(P) P$x)),
            max(sapply(posterior[id], function(P) P$x)))
  plot(posterior[[id[1]]], xlab=parname, xlim=xlim,
       )
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





