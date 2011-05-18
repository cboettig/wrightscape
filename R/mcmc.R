#multiou can try and take lca as a parameter option rather than calculating each time, for efficiency

step_fn <- function(pars, stepsizes = .02){
# Sequential random updating 
  j <- sample(1:length(pars), 1)
  pars[j] <- rnorm(1, pars[j], stepsizes)
  pars
}

## Proposal density is symmetric, so we won't need Q
Q <- function(pars, proposed){
  dnorm(pars, proposed, stepsizes=1, log=TRUE)
#    alpha <- loglik(pars) + prior(pars) + Q(pars, proposed) - loglik(proposed) - prior(proposed) - Q(proposed, pars)
}

mcmc_fn <- function(pars, loglik, prior, MaxTime=1e3, stepsizes=.02, ...){
  history <- matrix(NA, nrow=MaxTime, ncol=(1+length(pars)))
  for(t in 1:MaxTime){
    Pi <- loglik(pars)+prior(pars)
    history[t,] <- c(Pi, pars)
    proposed <- step_fn(pars, stepsizes)
    # Hastings ratio
    if (loglik(proposed)+prior(proposed) - Pi > runif(1) )
      pars <- proposed
  }
  history
}

beta <- function(i, Delta_T=1){
  1/(1+Delta_T*(i-1))
}


mcmcmc_fn <- function(pars, loglik, prior, MaxTime=1e3, indep=100, stepsizes=.02, ...){
# Metropolis Coupled Markov Chain Monte Carlo
# Args:
#   pars: a list of length n_chains, with elements that pars[[i]] that can be passed to loglik
#   loglik: a function to calculate the log-likelihood of chain i at pars[[i]], 
#   prior: a function to calculate the prior density
#   MaxTime: length of time to run the chain
  n_chains <- length(pars)
  n_pars <- length(pars[[1]])

# in case we want to store the complete history.  Should have the option
# of writing this to a file for speed? Or do that all in C...
  chains <- lapply(1:n_chains, function(i)  matrix(NA, nrow=MaxTime, ncol=(1+n_pars)) )

  Interval <- matrix(1:MaxTime, nrow=indep)

  for(s in 1:(MaxTime/indep)){
  # Evolve chains independently for "indep" time steps
    out <- sfLapply(1:n_chains, 
            function(i){
              out <- matrix(NA, ncol=n_pars+1, nrow=indep)
              for(t in 1:indep){
                Pi <- loglik(pars[[i]]) + prior(pars[[i]])
                out[t,] <- c(Pi, pars[[i]]) # more simply could print this step
                proposed <- step_fn(pars[[i]], stepsizes)
                if ( beta(i) * ( loglik(proposed)+prior(proposed) - Pi ) > runif(1) )
                  pars[[i]] <- proposed
              }
              out
            })
    # write to output
    for(i in 1:n_chains){
      chains[[i]][Interval[,s],] <- out[[i]]
      pars[[i]] <- out[[i]][indep,][-1] # copy parameters (but not loglik)
    }
    # Propose a chain swap every "indep" time steps beween pars_i and j
    pick <- sample(1:length(pars), 2)
    i <- pick[1]; j <- pick[2] # for convience
    R <- 
      beta(i) * ( loglik(pars[[j]]) - loglik(pars[[i]]) ) +
      beta(j) * ( loglik(pars[[i]]) - loglik(pars[[j]]) )

print(paste("swap chain", i, "with", j, "proposed, R =", exp(R)))

    if(exp(R) > runif(1)){
print("swap accepted")
      pars[[i]] <- pars[[j]] # accept the swap
    }

  }
  chains 
}


# 
rc_prior <- function(pars){
# Prior for global theta, global sigma, regime-specific alpha
  n_regimes <- length(pars) - 3  

  out <- dnorm(pars[1], 0, 1000, log=TRUE) +  # Xo prior -- 
  sum(dexp(pars[2:(1+n_regimes)], 1, log=TRUE)) +  # alpha prior -- 
  dexp(pars[2+n_regimes], 1, log=TRUE) + # sigma prior -- 
  dnorm(pars[3+n_regimes], 0, 1000, log=TRUE) # theta prior 
  out
}


ws_mcmc <- function(data, tree, regimes, alpha=1,
                    sigma=1, theta=NULL, Xo=NULL, prior=rc_prior, MaxTime = 1e3, ...){

## Initialize the parameters for the model
# alpha varies by regime, theta and sigma are global
# par is Xo, all alphas, theta, sigma
  n_regimes <- length(levels(regimes))
  par <- numeric(n_regimes+3)

  if(is.null(Xo))
    Xo <- mean(data, na.rm=TRUE) 
  if(is.null(theta))
    theta <- Xo 
  par[1] <- Xo
  if(length(alpha) == n_regimes){
      par[2:(1+n_regimes)] <- alpha
  } else {
      par[2:(1+n_regimes)] <- rep(alpha, n_regimes)
  }
  par[2+n_regimes] <- sigma 
  par[3+n_regimes] <- theta

  lca <- lca_calc(tree)

  # Likelihood as a function of optimizable parameters
  f <- function(par){
      Xo <- par[1]
      alpha <- par[2:(1+n_regimes)]
      sigma <- rep(par[2+n_regimes], n_regimes) 
      theta <- rep(par[3+n_regimes], n_regimes) 
      if (any(alpha < 0)){
          llik <- -Inf
      }
      else if (any(sigma<0)){
          llik <- -Inf
      } else {
          llik<-multiOU_lik_lca(data, tree, regimes, alpha=alpha,
                                sigma=sigma, theta=theta, Xo=Xo, lca)
      }
      llik
  }
 # sfLapply(1:4, function(i){
  #  par <- abs(rnorm(length(par), par, sd=3))
    mcmc_fn(par,f, prior, MaxTime=MaxTime, ...) 
 # })
}



general_prior <- function(pars){
# Prior for regime-specific alpha, theta, sigma
  n_regimes <- (length(pars) - 1)/3

  out <- dnorm(pars[1], 0, 1000, log=TRUE) +  # Xo prior -- 
  sum(dexp(pars[2:(1+n_regimes)], 1, log=TRUE)) +  # alpha prior -- 
  sum(dexp(pars[(2+n_regimes):(1+2*n_regimes)], 1, log=TRUE)) + # sigma prior -- 
  sum(dnorm(pars[(2+2*n_regimes):(1+3*n_regimes)], 0, 1000, log=TRUE)) # theta prior 
  out
}


general_mcmc <- function(data, tree, regimes, alpha=1, sigma=1, Xo=NULL, prior=general_prior, MaxTime, indep=MaxTime/100, ...){
# all are regime dependent
    # intialize a parameter vector to optimize: 
    # Xo, alpha, sigma, and the n_regime thetas
    n_regimes <- length(levels(regimes))
    par <- numeric(1+3*n_regimes)

    # Some starting conditions
    if(is.null(Xo)) Xo <- mean(data, na.rm=TRUE) 
    par[1] <- Xo
    if(length(alpha) == n_regimes){
        par[2:(1+n_regimes)] <- alpha
    } else {
        par[2:(1+n_regimes)] <- rep(alpha, n_regimes)
    }
    if(length(sigma) == n_regimes){
        par[(2+n_regimes):(1+2*n_regimes)] <- sigma 
    } else {
        par[(2+n_regimes):(1+2*n_regimes)] <- rep(sigma, n_regimes)
    }
    par[(2+2*n_regimes):(1+3*n_regimes)] <- rep(Xo, n_regimes)
    lca <- lca_calc(tree)
    # Likelihood as a function of optimizable parameters
    f <- function(par){
        Xo <- par[1]
        alpha <- par[2:(1+n_regimes)]
        sigma <- par[(2+n_regimes):(1+2*n_regimes)]
        theta <- par[(2+2*n_regimes):(1+3*n_regimes)]
        if (any(alpha < 0)){
            llik <- -Inf
        }
        else if (any(sigma<0)){
            llik <- -Inf
        } else {
            llik<-multiOU_lik_lca(data, tree, regimes, alpha=alpha,
                                  sigma=sigma, theta=theta, Xo=Xo, lca)
        }
        llik
    }
    par2 <- par; par3 <- par
    par2[2:(1+n_regimes)] = .001 # start with small alphas
    par3[2:(1+n_regimes)] = 5 # start with large alphas
    par4 <- sapply(rexp(length(par), 5), abs)


    pars <- list(par, par2, par3, par4)
    mcmcmc_fn(pars,f, prior, MaxTime=MaxTime, ...) 
}


