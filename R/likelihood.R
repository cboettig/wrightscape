#' Fits the generic multitype OU.  
#'
#' @details 
#' Submodels such as brownie, and other unique models, can be created
#' by specifiying how parameters are treated under model_spec

#' get the likelihood of the specified model using the specified parameters
#' @param data the trait data
#' @param tree the phylogenetic tree in ouch format
#' @param regimes the regimes in ouch format
#' @param model_spec a list that specifies the model, see details
#' @param Xo root state
#' @param alpha a vector of length n_regimes if indep in model, or a scalar
#' @param sigma a vector of length n_regimes if indep in model, or a scalar
#' @param theta  a vector of length n_regimes if indep in model, or a scalar
#' @return log likelihood
#' @details
#' the general model form is specified by model_spec list.  This specifies 
#' which parameters out of alpha, theta, and sigma are independently estimated
#' on each regime, kept global across regimes, or, in the case of alpha,
#' fixed to zero (to give purely Brownian behavior).  i.e.
#' ouch model is equivalent to: list(alpha="global", sigma="global",
#' theta="indep"), while the brownie model is equivalent to 
#' list(alpha="fixed", sigma="indep", theta="global") 
#' @useDynLib wrightscape
#' @import ouch
#' @import geiger
#' @export
multiTypeOU <- function(data, tree, regimes, model_spec =
                       list(alpha="indep", sigma="indep", theta="indep"),
                       Xo=NULL, alpha=1, sigma=1, theta=NULL, ...){
  opts=list(...)

  myCall <- match.call()
  n_regimes <- length(levels(regimes))
  par <- setup_pars(data, tree, regimes, model_spec, Xo=Xo, 
                    alpha=alpha, sigma=sigma, theta=theta)
  f <- llik.closure(data, tree, regimes, model_spec, neg=TRUE)
  optim_output <- optim(par,f, ...) 
  indices <- get_indices(model_spec, n_regimes)

  ## deal with fixed alpha at 0
  if(is.null(indices$alpha_i))
    alpha_out <- rep(1e-12, n_regimes)
  else 
    alpha_out <- optim_output$par[indices$alpha_i]
  Xo_out <- optim_output$par[indices$theta_i][match(regimes[1], levels(regimes))]
  output <- list(data=data, tree=tree, regimes=regimes, 
                 loglik=-optim_output$value, 
                 alpha=alpha_out, 
                 sigma=optim_output$par[indices$sigma_i],
                 theta=optim_output$par[indices$theta_i],
                 Xo = Xo_out,
                 optim_output=optim_output, model_spec=model_spec,
                 convergence=optim_output$convergence,
                 opts=opts)
  class(output) = "multiOU"
  output

}

#' get the likelihood of the specified model using the specified parameters
#' @param data the trait data
#' @param tree the phylogenetic tree in ouch format
#' @param regimes the regimes in ouch format
#' @param model_spec a list that specifies the model, see details
#' @param Xo root state
#' @param alpha a vector of length n_regimes if indep in model, or a scalar
#' @param sigma a vector of length n_regimes if indep in model, or a scalar
#' @param theta  a vector of length n_regimes if indep in model, or a scalar
#' @return log likelihood
#' @details
#' the general model form is specified by model_spec list.  This specifies 
#' which parameters out of alpha, theta, and sigma are independently estimated
#' on each regime, kept global across regimes, or, in the case of alpha,
#' fixed to zero (to give purely Brownian behavior).  i.e.
#' ouch model is equivalent to: list(alpha="global", sigma="global",
#' theta="indep"), while the brownie model is equivalent to 
#' list(alpha="fixed", sigma="indep", theta="global") 
#' @export
evaluate_likelihood <- function(data, tree, regimes, model_spec =
                       list(alpha="indep", sigma="indep", theta="indep"),
                       Xo, alpha, sigma, theta){

  n_regimes <- length(levels(regimes))
  par <- setup_pars(data, tree, regimes, model_spec, Xo=Xo, 
                    alpha=alpha, sigma=sigma, theta=theta)
  f <- llik.closure(data, tree, regimes, model_spec, neg=FALSE)
  f(par)

}


#' returns the likelihood function for the chosen submodel
#' @keywords internal
llik.closure <- function(data, tree, regimes, model_spec, fixed=
                         list(alpha=1e-12), neg_llik=FALSE){
  n_regimes <- length(levels(regimes))
  indices <- get_indices(model_spec, n_regimes)
  lca <- lca_calc(tree)
  f <- function(par){
     Xo <- par[indices$theta_i][match(regimes[1], levels(regimes))]  # assumes root theta
    if(any(is.null(indices$alpha_i)))
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
    if(is.nan(llik))
      llik <- -Inf
    if(neg_llik) # for optimizer routines
      llik <- -llik
    llik
  }
  f
}


#' Helper function to work with llik.closure indexing 
#' @keywords internal
get_indices <- function(model_spec, n_regimes){
# Get the indices of the parameter vector (thing passed to the likelihood routine)
# that correspond to the different parameters.  Inverse of setup_pars function.  
# Examples
#   ouch: get_indices(list(alpha="global", sigma="global", theta="indep"), 2)
  alpha_start <- 1 # Xo is not to be estimated independently. 2 otherwise
  alpha_i <-  switch(model_spec$alpha, 
                     indep = alpha_start:(alpha_start+n_regimes-1), 
                     global = rep(alpha_start, n_regimes),
                     fixed = NULL)
  if(any(is.null(alpha_i))) 
    sigma_start <- alpha_start
  else 
    sigma_start <- max(alpha_i)+1
  sigma_i <-  switch(model_spec$sigma, 
                     indep = sigma_start:(sigma_start+n_regimes-1),
                     global = rep(sigma_start, n_regimes))
  if(any(is.null(sigma_i)))
    theta_start <- sigma_start
  else
    theta_start <- max(sigma_i)+1
  theta_i <-  switch(model_spec$theta, 
                     indep = theta_start:(theta_start+n_regimes-1), 
                     global = rep(theta_start,n_regimes))
  n <- max(theta_i) # index of last parameter is number of parameters
list(alpha_i=alpha_i, sigma_i=sigma_i, theta_i=theta_i, n=n)
}

#' Helper function to work with llik.closure indexing 
#' @keywords internal
setup_pars <- function(data, tree, regimes, model_spec, Xo=NULL, alpha=1, 
                       sigma=1, theta=NULL){
## Create the parameter vector matching model_spec from specified alpha, theta, sigma
## This is essentially the inverse function of get_indices
  n_regimes <- length(levels(regimes))
  indices <- get_indices(model_spec, n_regimes)
  pars <- numeric(indices$n)
  if(is.null(theta)) 
    theta <- mean(data, na.rm=TRUE)
  # pars[1] <- Xo # not to be estimated!
  if(!any(is.null(indices$alpha_i)))
    pars[indices$alpha_i] <- alpha
  pars[indices$sigma_i] <- sigma
  pars[indices$theta_i] <- theta 
  pars
}


#' Wrappers for the C functions that actually calculate the likelihood
#' @keywords internal
lca_calc <- function(tree){
# Calculates the last common ancestor matrix, which is used by the likelihood algorithm
# Because this is constant given a tree, it need not be recalculated each time.  
# This uses a rather naive C function 
# Args: 
#   tree: ouch-formatted tree
# Returns: 
#   A integer array (matrix) which can be passed to the likelihood function
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

#' compute the likelihood by passing data to the C level function, wrightscape 
#' @param data - ouch-style data
#' @param tree - ouch-tree
#' @param regimes - painting of selective regimes, as in ouch
#' @param alpha - (vector length n_regimes) gives strength of selection in each regime
#' @param sigma - (vector length n_regimes) gives diversification rate in each regime
#' @param theta - (vector length n_regimes) gives optimum trait in each regime
#' @param Xo - root value
#' @param lca - least common ancestor matrix, from lca_calc fun
#' @return the log likelihood at the given parameter values
#' @details 
#' A version of the multitype OU likelihood function that accepts the least 
#' common ancestor matrix as a parameter.  This may increase computational speed,
#' since the calculation needs to be done only once and is rather slow as implemented    
#' .C calls should do some error checking on the length of inputs maybe, 
#' to avoid crashes when given inappropriate calls
multiOU_lik_lca <- function(data, tree, regimes, alpha=NULL, sigma=NULL, theta=NULL, Xo=NULL, lca){

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













# A function that tries to guess smart starting conditions, probably overkill
smart_multiType <- function(data, tree, regimes, model_spec =
                       list(alpha="indep", sigma="indep", theta="indep"),
                       Xo=NULL, alpha=1, sigma=1, theta=NULL, ...){
## Fits the general multitype OU by seeding starting conditions from submodels
## Update will repeat all this crap, so go easy on reps if using SANN 
  opts=list(...)
  myCall <- match.call()
  n_regimes <- length(levels(regimes))

  ## Should the guesses be propigated through?  currently not. 
  alpha_spec = list(alpha="indep",  sigma="global", theta="global")   
  par <- setup_pars(data, tree, regimes, alpha_spec, Xo=Xo, 
                    alpha=alpha, sigma=sigma, theta=theta)
  alpha_f <- llik.closure(data, tree, regimes, alpha_spec, neg=TRUE)
  alpha_optim <- optim(par,alpha_f, ...) 
  alpha_indices <- get_indices(alpha_spec, n_regimes)

  sigma_spec = list(alpha="global", sigma="indep",  theta="global")   
  par <- setup_pars(data, tree, regimes, sigma_spec, Xo=Xo, 
                    alpha=alpha, sigma=sigma, theta=theta)
  sigma_f <- llik.closure(data, tree, regimes, sigma_spec, neg=TRUE)
  sigma_optim <- optim(par,sigma_f, ...) 
  sigma_indices <- get_indices(sigma_spec, n_regimes)

  theta_spec = list(alpha="global", sigma="global", theta="indep")   
  par <- setup_pars(data, tree, regimes, theta_spec, Xo=Xo, 
                    alpha=alpha, sigma=sigma, theta=theta)
  theta_f <- llik.closure(data, tree, regimes, theta_spec, neg=TRUE)
  theta_optim <- optim(par,theta_f, ...) 
  theta_indices <- get_indices(theta_spec, n_regimes)


  # get alpha guesses as the optimum alphas in the alpha model, etc
  alpha <- alpha_optim$par[alpha_indices$alpha_i]
  sigma <- sigma_optim$par[sigma_indices$sigma_i]
  theta <- theta_optim$par[theta_indices$theta_i]


  print(paste("alpha: ", alpha))
  print(paste("sigma: ", sigma))
  print(paste("theta: ", theta))
  print(paste("alpha LL", -alpha_optim$value, 
              "simga LL", -sigma_optim$value,
              "theta LL", -theta_optim$value))

  par <- setup_pars(data, tree, regimes, model_spec, Xo=Xo, 
                    alpha=alpha, sigma=sigma, theta=theta)
  f <- llik.closure(data, tree, regimes, model_spec, neg=TRUE)
  optim_output <- optim(par,f, ...) 
  indices <- get_indices(model_spec, n_regimes)


  ## deal with fixed alpha at 0
  if(is.null(indices$alpha_i))
    alpha_out <- rep(1e-12, n_regimes)
  else 
    alpha_out <- optim_output$par[indices$alpha_i]

  output <- list(data=data, tree=tree, regimes=regimes, 
                 loglik=-optim_output$value, Xo=optim_output$par[1], 
                 alpha = alpha_out, 
                 sigma=optim_output$par[indices$sigma_i],
                 theta=optim_output$par[indices$theta_i],
                 optim_output=optim_output, model_spec=model_spec,
                 convergence=optim_output$convergence,
                 opts=opts)
  class(output) = "multiOU"
  output

}


