#multiou can try and take lca as a parameter option rather than calculating each time, for efficiency

step_fn <- function(pars, stepsizes = .1){
  rnorm(1, sample(pars, 1), stepsizes)
  #rnorm(length(pars), pars, stepsizes)
}

Q <- function(pars, proposed){
  dnorm(pars, proposed, stepsizes=1, log=TRUE)
}
#    alpha <- loglik(pars) + Q(pars, proposed) - loglik(proposed) - Q(proposed, pars)

mcmc <- function(pars, loglik, MaxTime=1e3, ...){

  history <- matrix(NA, nrow=MaxTime, ncol=(1+length(pars)))
  for(t in 1:MaxTime){
    L <- loglik(pars)
    history[t,] <- c(L, pars)
    proposed <- step_fn(pars)
    # Hastings ratio
    if (loglik(proposed) - L > runif(1) )
      pars <- proposed
  }

  history
}

ws_mcmc <- function(data, tree, regimes, alpha=NULL,
                    sigma=NULL, theta=NULL, Xo=NULL, MaxTime = 1e3, ...){

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


  mcmc(par,f, MaxTime=MaxTime, ...) 
}


