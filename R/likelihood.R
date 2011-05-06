#multiou can try and take lca as a parameter option rather than calculating each time, for efficiency

update.multiOU <- function(model, data){
    switch(model$submodel,
           wright = wright(data, model$tree, model$regimes,
                           model$Xo, model$alpha, model$sigma),
           ouch = ouch(data, model$tree, model$regimes,
                       Xo=model$Xo, alpha=model$alpha, sigma=model$sigma),
           brownie = brownie(data, model$tree, model$regimes,
                             model$sigma))
}

release_constraint <- function(data, tree, regimes, alpha=NULL,
                               sigma=NULL, theta=NULL, Xo=NULL, ...){
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
          llik <- -1e12
      }
      else if (any(sigma<0)){
          llik <- -1e12
      } else {
          llik<-multiOU_lik_lca(data, tree, regimes, alpha=alpha,
                                sigma=sigma, theta=theta, Xo=Xo, lca)
      }
      -llik
  }
  optim_output <- optim(par,f, ...) 
#    optim(par,f, method="L", lower=c(-Inf, rep(0,n_regimes), rep(-Inf, n_regimes), rep(0, n_regimes))) 
  output <- list(data=data, tree=tree, regimes=regimes, 
                 loglik=-optim_output$value, Xo=optim_output$par[1], 
                 alpha=optim_output$par[2:(1+n_regimes)], 
                 theta=optim_output$par[2+n_regimes],
                 sigma=optim_output$par[3+n_regimes],
                 optim_output=optim_output, submodel="release_constraint",
                 convergence=optim_output$convergence)
  class(output) = "multiOU"
  output
}


wright <- function(data, tree, regimes, alpha=1, sigma=1, Xo=NULL, ...){
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
            llik <- -1e12
        }
        else if (any(sigma<0)){
            llik <- -1e12
        } else {
            llik<-multiOU_lik_lca(data, tree, regimes, alpha=alpha,
                                  sigma=sigma, theta=theta, Xo=Xo, lca)
        }
        -llik
    }
    optim_output <- optim(par,f, ...) 
#    optim(par,f, method="L", lower=c(-Inf, rep(0,n_regimes), rep(-Inf, n_regimes), rep(0, n_regimes))) 
    output <- list(data=data, tree=tree, regimes=regimes, 
                   loglik=-optim_output$value, Xo=optim_output$par[1], 
                   alpha=optim_output$par[2:(1+n_regimes)], 
                   theta=optim_output$par[(2+2*n_regimes):(1+3*n_regimes)],
                   sigma=optim_output$par[(2+n_regimes):(1+2*n_regimes)],
                   optim_output=optim_output, submodel="wright",
                   convergence=optim_output$convergence)
    class(output) = "multiOU"
    output
}




# OUCH
ouch <- function(data, tree, regimes, alpha=1, sigma=1, Xo=NULL, ...){
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
                   convergence=optim_output$convergence)
    class(output) = "multiOU"
    output
}

# Brownie
brownie <- function(data, tree, regimes, sigma=1, ...){ 

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
        alpha <- rep(1e-12, n_regimes) 
        theta <- rep(Xo, n_regimes)
        llik <- multiOU_lik_lca(data, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo, lca)
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
                   convergence=optim_output$convergence)
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
