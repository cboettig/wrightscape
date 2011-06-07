

#multiou can try and take lca as a parameter option rather than calculating each time, for efficiency
update.multiOU <- function(model, data){
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


######## Use these to define a generic model ######## 
get_indices <- function(model_spec, n_regimes){
# Examples
#   ouch: get_indices(list(alpha="global", sigma="global", theta="indep"), 2)
  alpha_start <- 2
  alpha_i <-  switch(model_spec$alpha, 
                     indep = alpha_start:(alpha_start+n_regimes-1), 
                     global = rep(alpha_start, n_regimes),
                     fixed = NA)
  if(any(is.na(alpha_i))) 
    sigma_start <- alpha_start
  else 
    sigma_start <- max(alpha_i)+1
  sigma_i <-  switch(model_spec$sigma, 
                     indep = sigma_start:(sigma_start+n_regimes-1),
                     global = rep(sigma_start, n_regimes))
  if(any(is.na(sigma_i)))
    theta_start <- sigma_start
  else
    theta_start <- max(sigma_i)+1
  theta_i <-  switch(model_spec$theta, 
                     indep = theta_start:(theta_start+n_regimes-1), 
                     global = rep(theta_start,n_regimes))
  n <- max(theta_i)
list(alpha_i=alpha_i, sigma_i=sigma_i, theta_i=theta_i, n=n)
}


setup_pars <- function(data, tree, regimes, model_spec, Xo=NULL, alpha=1, sigma=1, theta=NULL){
  n_regimes <- length(levels(regimes))
  indices <- get_indices(model_spec, n_regimes)
  pars <- numeric(indices$n)
  if(is.null(Xo)) 
    Xo <- mean(data, na.rm=TRUE) # should use phylo mean, but need to convert tree to ape type
  if(is.null(theta)) 
    theta <- mean(data, na.rm=TRUE)
  pars[1] <- Xo
  if(!any(is.na(indices$alpha_i)))
    pars[indices$alpha_i] <- alpha
  pars[indices$sigma_i] <- sigma
  pars[indices$theta_i] <- theta 
  pars
}



llik.closure <- function(data, tree, regimes, model_spec, fixed=list(alpha=1e-12), neg_llik=FALSE){
  n_regimes <- length(levels(regimes))
  indices <- get_indices(model_spec, n_regimes)
  lca <- lca_calc(tree)
  f <- function(par){
    Xo <- par[1]
    if(any(is.na(indices$alpha_i)))
      alpha <- rep(fixed$alpha, n_regimes)
    else 
      alpha <- par[indices$alpha_i]
    sigma <- par[indices$sigma_i]
    theta <- par[indices$theta_i] 
    if (any(alpha < 0))
        llik <- -Inf
    else if (any(sigma<0))
        llik <- -Inf
    else 
        llik<-multiOU_lik_lca(data, tree, regimes, alpha=alpha,
                              sigma=sigma, theta=theta, Xo=Xo, lca)
    if(neg_llik) # for optimizer routines
      llik <- -llik
    llik
  }
  f
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


lca_calc <- function(tree){
# Calculates the last common ancestor matrix, which is used by the likelihood algorithm
# Because this is constant given a tree, it need not be recalculated each time.  
# This uses a rather naive C function 
#
# Args: tree ouch-formatted tree
#
# Returns: A integer array (matrix) which can be passed to the likelihood function
#

#  Format and error handling -- FIXME: this should be generalized for these functions

    if (is(tree, "phlyo")) tree <- ape2ouch(tree)
    if (!is(tree, "ouchtree")) warning("tree is not a ouch-tree format")
    ## CONVERT tree elements into C formats 
	ancestor <- as.numeric(tree@ancestors)
	ancestor[is.na(ancestor)] = 0 
	ancestor <- ancestor-1  # C-style indexing
	## ouch gives cumulative time, not branch-length!!
	anc <- as.integer(tree@ancestors[!is.na(tree@ancestors)])
	lengths <- c(0, tree@times[!is.na(tree@ancestors)] - tree@times[anc] )
	branch_length <- lengths/max(tree@times)
	n_nodes <- length(branch_length)
    ## calculate the lca_matrix
    lca_matrix <- integer(n_nodes^2)
    lca <- .C("calc_lca", as.integer(ancestor), as.double(branch_length),
              as.integer(n_nodes), as.integer(lca_matrix))
    lca[[4]]
}

## .C calls should do some error checking on the length of inputs maybe, to avoid crashes when given inappropriate calls
multiOU_lik_lca <- function(data, tree, regimes, alpha=NULL, sigma=NULL, theta=NULL, Xo=NULL, lca){
# A version of the multitype OU likelihood function that accepts the least 
# common ancestor matrix as a parameter.  This may increase computational speed,
# since the calculation needs to be done only once and is rather slow as implemented    
# 
# Args: data -- ouch-style data
#       tree -- ouch-tree
#       regimes -- painting of selective regimes, as in ouch
#       alpha -- (vector length n_regimes) gives strength of selection in each regime
#       sigma -- (vector length n_regimes) gives diversification rate in each regime
#       theta -- (vector length n_regimes) gives optimum trait in each regime
#       Xo -- root value
#       lca -- least common ancestor matrix, from lca_calc fun
# Returns: log likelihood 

    ## ERROR HANDLING, write this as a seperate function
	# data should be a numeric instead of data.frame.  Should check node names or node order too!
	dataIn <- data
	if(is(data, "data.frame") | is(data, "list")) { 
		data <- data[[1]]
		if( !is(data, "numeric")) {stop("data should be data frame or numeric") }
	}
	# regimes should be a factor instead of data.frame
	regimesIn <- regimes
	if(is(regimes, "data.frame")) { 
		regimes <- regimes[[1]]
		if( !is(regimes, "factor")) {stop("unable to interpret regimes") }
	}

    ## CONVERT tree elements into C formats 
	if(is.null(Xo)){ Xo <- mean(data, na.rm=TRUE) }
	data[is.na(data)] = 0 
	ancestor <- as.numeric(tree@ancestors)
	ancestor[is.na(ancestor)] = 0 
	ancestor <- ancestor-1  # C-style indexing
	## ouch gives cumulative time, not branch-length!!
	anc <- as.integer(tree@ancestors[!is.na(tree@ancestors)])
	lengths <- c(0, tree@times[!is.na(tree@ancestors)] - tree@times[anc] )
	branch_length <- lengths/max(tree@times)
	n_nodes <- length(branch_length)
	n_regimes <- length(levels(regimes))
	if(length(alpha) == 1){ alpha <- rep(alpha, n_regimes) }
	if(is.null(theta)) { theta <- rep(Xo, n_regimes) }
	if(length(sigma) == 1) { sigma <- rep(sigma, n_regimes) }
	levels(regimes) <- 1:n_regimes
	regimes <- as.integer(regimes)-1  # convert to C-style indexing

    llik <- 0
    ## calculate the likelihood
	o<- .C("calc_lik", as.double(Xo), as.double(alpha), as.double(theta),
            as.double(sigma), as.integer(regimes), as.integer(ancestor),
            as.double(branch_length), as.double(data), as.integer(n_nodes),
            as.integer(lca), as.double(llik))

# loglikelhood
	o[[11]] 
}








## Recalculates lca each time, useful for a single call, 
# but too slow to use in optimization.  
# Currently exported, while separate lca and multiOU_lik functions are not.  
multiOU_lik <- function(data, tree, regimes, alpha=NULL, sigma=NULL, theta=NULL, Xo=NULL){

    ## ERROR HANDLING, write this as a seperate function
	# data should be a numeric instead of data.frame.  Should check node names or node order too!
	dataIn <- data
	if(is(data, "data.frame") | is(data, "list")) { 
		data <- data[[1]]
		if( !is(data, "numeric")) {stop("data should be data frame or numeric") }
	}
	# regimes should be a factor instead of data.frame
	regimesIn <- regimes
	if(is(regimes, "data.frame")) { 
		regimes <- regimes[[1]]
		if( !is(regimes, "factor")) {stop("unable to interpret regimes") }
	}

    ## CONVERT tree elements into C formats 
	if(is.null(Xo)){ Xo <- mean(data, na.rm=TRUE) }
	data[is.na(data)] = 0 
	ancestor <- as.numeric(tree@ancestors)
	ancestor[is.na(ancestor)] = 0 
	ancestor <- ancestor-1  # C-style indexing
	## ouch gives cumulative time, not branch-length!!
	anc <- as.integer(tree@ancestors[!is.na(tree@ancestors)])
	lengths <- c(0, tree@times[!is.na(tree@ancestors)] - tree@times[anc] )
	branch_length <- lengths/max(tree@times)
	n_nodes <- length(branch_length)
	n_regimes <- length(levels(regimes))
	if(length(alpha) == 1){ alpha <- rep(alpha, n_regimes) }
	if(is.null(theta)) { theta <- rep(Xo, n_regimes) }
	if(length(sigma) == 1) { sigma <- rep(sigma, n_regimes) }

	levels(regimes) <- 1:n_regimes
	regimes <- as.integer(regimes)-1  # convert to C-style indexing



    ## calculate the lca_matrix
    lca_matrix <- integer(n_nodes^2)
    lca <- .C("calc_lca", as.integer(ancestor), as.double(branch_length),
              as.integer(n_nodes), as.integer(lca_matrix))

    llik <- 0
    ## calculate the likelihood
	o<- .C("calc_lik", as.double(Xo), as.double(alpha), as.double(theta),
            as.double(sigma), as.integer(regimes), as.integer(ancestor),
            as.double(branch_length), as.double(data), as.integer(n_nodes),
            as.integer(lca[[4]]), as.double(llik))

	output <- list(data=dataIn, tree=tree, regimes=regimesIn, loglik=o[[11]], Xo=o[[1]],
                   alpha = o[[2]], theta =  o[[3]], sigma = o[[4]]  )  
	class(output) <- "multiOU"
	output
}




