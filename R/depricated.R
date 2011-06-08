# depricated functions that still have dependencies in the demos
# brownie, OUCH, wright, release_constraint 


#multiou can try and take lca as a parameter option rather than calculating each time, for efficiency
old_update.multiOU <- function(model, data){
    switch(model$submodel,
           wright = do.call(wright, 
                            c(list(data=data, tree=model$tree,
                                   regimes=model$regimes,
                                   Xo=model$Xo, alpha=model$alpha, 
                                  sigma=model$sigma, theta=model$theta),
                            model$opts)),
           ouch = do.call(ouch, 
                          c(list(data=data, tree=model$tree, 
                                 regimes=model$regimes, Xo=model$Xo, 
                                 alpha=model$alpha, sigma=model$sigma),
                          model$opts)),
           brownie = do.call(brownie, 
                             c(list(data=data, tree=model$tree,
                                    regimes=model$regimes, sigma=model$sigma),
                             model$opts)),
           release_constraint =  
            do.call(release_constraint, 
                    c(list(data=data, tree=model$tree,
                           regimes=model$regimes, Xo=model$Xo, 
                           alpha=model$alpha, sigma=model$sigma,
                           theta=model$theta), 
                    model$opts)))
           
}


# Likelihood as a function of optimizable parameters
llik.release = function(data, tree, regimes, lca=NULL){
# returns a likelihood funciton of pars: {Xo, alpha, sigma, theta}
  n_regimes <- length(levels(regimes))
  if(is.null(lca))
    lca <- lca_calc(tree)
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
      -llik
  }
  f
}


release_constraint <- function(data, tree, regimes, alpha=NULL,
                               sigma=NULL, theta=NULL, Xo=NULL, ...){
  opts <- list(...)
# alpha varies by regime, theta and sigma are global
# par is Xo, all alphas, theta, sigma

## Probably a much better way to handle this.  could also be functionalized...

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

  f <- llik.release(data, tree, regimes, lca)
  print(par)
  print(paste("starting loglik = ", -f(par)))


  optim_output <- optim(par,f, ...) 
#    optim(par,f, method="L", lower=c(-Inf, rep(0,n_regimes), rep(-Inf, n_regimes), rep(0, n_regimes))) 
  output <- list(data=data, tree=tree, regimes=regimes, 
                 loglik=-optim_output$value, Xo=optim_output$par[1], 
                 alpha=optim_output$par[2:(1+n_regimes)], 
                 sigma=optim_output$par[2+n_regimes],
                 theta=optim_output$par[3+n_regimes],
                 optim_output=optim_output, submodel="release_constraint",
                 convergence=optim_output$convergence,
                 opts=opts)
  class(output) = "multiOU"
  output
}



# Likelihood as a function of optimizable parameters
llik.wright <- function(data, tree, regimes, lca=NULL){
  n_regimes <- length(levels(regimes))

  if(is.null(lca))
    lca <- lca_calc(tree)
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
    -llik
  }
  f
}


wright <- function(data, tree, regimes, alpha=1, sigma=1, Xo=NULL, theta=NULL, ...){
  opts <- list(...)


    # all are regime dependent
    # intialize a parameter vector to optimize: 
    # par = {Xo, alphas, sigmas, thetas}
    n_regimes <- length(levels(regimes))
    par <- numeric(1+3*n_regimes)

    # Create par vector from the given starting conditions 
    # (or guess if not given)

    ## Specify the indices of each. Xo is always 1
    alpha_i <- 2:(1+n_regimes) 
    sigma_i <- (2+n_regimes):(1+2*n_regimes)
    theta_i <- (2+2*n_regimes):(1+3*n_regimes)

    if(is.null(Xo)) Xo <- mean(data, na.rm=TRUE) 
    par[1] <- Xo
    if(length(alpha) == n_regimes){
        par[alpha_i] <- alpha
    } else {
        par[alpha_i] <- rep(alpha, n_regimes)
    }
    if(length(sigma) == n_regimes){
        par[sigma_i] <- sigma 
    } else {
        par[sigma_i] <- rep(sigma, n_regimes)
    } 
    if(is.null(theta)){
      par[theta_i] <- rep(Xo, n_regimes)
    } else if(length(theta) == n_regimes){
       par[theta_i] <- theta
    } else {
      par[theta_i] <- rep(theta, n_regimes)
    }
    lca <- lca_calc(tree)
    # Likelihood as a function of optimizable parameters
    f <- llik.wright(data, tree, regimes, lca)

    print(par)
    print(paste("starting loglik = ", -f(par)))


    optim_output <- optim(par,f, ...) 
#    optim(par,f, method="L", lower=c(-Inf, rep(0,n_regimes), rep(-Inf, n_regimes), rep(0, n_regimes))) 
    output <- list(data=data, tree=tree, regimes=regimes, 
                   loglik=-optim_output$value, Xo=optim_output$par[1], 
                   alpha=optim_output$par[alpha_i], 
                   sigma=optim_output$par[sigma_i],
                   theta=optim_output$par[theta_i],
                   optim_output=optim_output, submodel="wright",
                   convergence=optim_output$convergence, opts=opts)
    class(output) = "multiOU"
    output
}


# OUCH
ouch <- function(data, tree, regimes, alpha=1, sigma=1, Xo=NULL, ...){
  opts <- list(...)

# alpha is fixed at ~zero, sigma is regime dependent, theta is global

    # intialize a parameter vector to optimize: 
    # Xo, alpha, sigma, and the n_regime thetas
    n_regimes <- length(levels(regimes))
    par <- numeric(3+n_regimes)

    if(length(alpha) > 1){
      alpha <- alpha[1]
    }
    if(length(sigma) > 1){
      sigma <- sigma[1]
    }

    # Some starting conditions
    if(is.null(Xo)) Xo <- mean(data, na.rm=TRUE) 
    par[1] <- Xo
    par[2] <- alpha
    par[3] <- sigma
    par[4:(3+n_regimes)] <- rep(Xo, n_regimes)
    lca <- lca_calc(tree)

    # Likelihood as a function of optimizable parameters
    f <- function(par){
        Xo <- par[1]
        alpha <- rep(par[2], n_regimes)
        theta <- par[4:(3+n_regimes)]
        sigma <- rep(par[3], n_regimes) # everything else
        if (any(alpha < 0)){ 
            llik <- -Inf
        }
        else if (any(sigma<0)){
            llik <- -Inf
        } else {
            llik <- multiOU_lik_lca(data, tree, regimes, alpha=alpha,
                                    sigma=sigma, theta=theta, Xo=Xo, lca)
        }
        -llik
    }
    optim_output <- optim(par,f, ...) 
    output <- list(data=data, tree=tree, regimes=regimes, 
                   loglik=-optim_output$value, Xo=optim_output$par[1], 
                   alpha=optim_output$par[2], 
                   theta=optim_output$par[4:(3+n_regimes)],
                   sigma=optim_output$par[3],
                   optim_output=optim_output,
                   submodel="ouch",
                   convergence=optim_output$convergence, opts=opts)
    class(output) = "multiOU"
    output
}


# Brownie
# should take Xo
brownie <- function(data, tree, regimes, sigma=1, ...){ 
  opts <- list(...)

    # intialize a parameter vector to optimize: 
    # Xo, followed by the n_regime sigmas
    n_regimes <- length(levels(regimes))
    pars <- numeric(1+n_regimes)
    # Some starting conditions
    pars[1] <- mean(data, na.rm=TRUE) #Xo
    if(length(sigma) == n_regimes){
        pars[2:(1+n_regimes)] <- sigma # sigmas
    } else {
        pars[2:(1+n_regimes)] <- rep(sigma, n_regimes) # sigmas
    }
    lca <- lca_calc(tree)

    # Likelihood as a function of optimizable parameters
    f <- function(pars){
        Xo <- pars[1]
        sigma <- pars[2:(1+n_regimes)] # everything else
        alpha <- rep(1e-12, n_regimes) ## all alphas approx 0
        theta <- rep(Xo, n_regimes)
        if (any(sigma<0)){
            llik <- -Inf
        } else {
        llik <- multiOU_lik_lca(data, tree, regimes, alpha=alpha, sigma=sigma,
                                theta=theta, Xo=Xo, lca)
        }
        -llik
    }
    optim_output <- optim(pars,f, ...) 
#    optim(par,f, method="L", lower=c(-Inf, rep(0,n_regimes), rep(-Inf, n_regimes), rep(0, n_regimes))) 
    output <- list(data=data, tree=tree, regimes=regimes, 
                   loglik=-optim_output$value, Xo=optim_output$par[1], 
                   alpha=rep(1e-12, n_regimes),
                   theta=rep(pars[1], n_regimes),
                   sigma=optim_output$par[2:(1+n_regimes)],
                   optim_output=optim_output,
                   submodel="brownie",
                   convergence=optim_output$convergence, opts=opts)
    class(output) = "multiOU"
    output
}


