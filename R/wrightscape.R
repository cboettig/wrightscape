# File: wrightscape.R
# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD

# simulation function, calls the C code
# @keywords internal 
simulate_wrightscape <- function(tree, regimes, Xo, alpha, theta, sigma,
                                 seed=NULL){

	# regimes should be a factor instead of data.frame
	regimesIn <- regimes
	if(is(regimes, "data.frame")) { 
		regimes <- regimes[[1]]
		if( !is(regimes, "factor")) {stop("unable to interpret regimes") }
	}

	ancestor <- as.numeric(tree@ancestors)
	ancestor[is.na(ancestor)] = 0 
	ancestor <- ancestor-1  # C-style indexing

	## ouch gives cumulative time, not branch-length!!
	anc <- as.integer(tree@ancestors[!is.na(tree@ancestors)])
	lengths <- c(0, tree@times[!is.na(tree@ancestors)] - tree@times[anc] )
	branch_length <- lengths/max(tree@times)

	n_nodes <- length(branch_length)
	n_regimes <- length(levels(regimes))
	n_tips <- (n_nodes+1)/2

	levels(regimes) <- 1:n_regimes
	regimes <- as.integer(regimes)-1  # convert to C-style indexing

  if(is.null(seed))
  	seed <- runif(1)*2^16

	if(length(alpha) == 1){ alpha <- rep(alpha, n_regimes) }
	if(is.null(theta)) { theta <- rep(Xo, n_regimes) }
	if(length(sigma) == 1) { sigma <- rep(sigma, n_regimes) }

## Should make sure return values of internal nodes are NAs (unless requested)
## should also make sure outputs aren't NANs!

	o<- .C("simulate_model",
		as.double(Xo),
		as.double(alpha),
		as.double(theta),
		as.double(sigma),
		as.integer(regimes),
		as.integer(ancestor),
		as.double(branch_length),
		double(n_nodes),
		as.integer(n_nodes),
		as.integer(n_regimes),
		double(1), 
		as.double(seed)
	  )
	simdata <- data.frame(o[[8]], row.names = tree@nodes)
	output <- list(rep.1=simdata, tree=tree, regimes=regimesIn, loglik=o[[11]],
                 alpha = o[[2]], theta =  o[[3]], sigma = o[[4]]  )  
	class(output) <- "wrighttree"
	output
}


#' fit the full model entirely in C code, for speed. 
#' @details
#' Does not allow fitting of the submodels.  Largely historical
#' @keywords internal 
wrightscape <- function(data, tree, regimes, alpha=1, sigma=1, 
                        theta = NULL, Xo = NULL, use_siman=0){

	# data should be a numeric instead of data.frame.  
  # Should check node names or node order too!
	dataIn <- data
	if(is(data, "data.frame") | is(data, "list")) { 
		data <- data[[1]]
		if( !is(data, "numeric")) {
      stop("data should be data frame or numeric") 
    }
	}

	# regimes should be a factor instead of data.frame
	regimesIn <- regimes
	if(is(regimes, "data.frame")) { 
		regimes <- regimes[[1]]
		if( !is(regimes, "factor")) {
      stop("unable to interpret regimes") 
    }
	}


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

  # rep these if not specified
	if(length(alpha) == 1){ alpha <- rep(alpha, n_regimes) }
	if(is.null(theta)) { theta <- rep(Xo, n_regimes) }
	if(length(sigma) == 1) { sigma <- rep(sigma, n_regimes) }

  ## shouldn't be necessary
  if(is.list(alpha)) unlist(alpha)
  if(is.list(theta)) unlist(theta)
  if(is.list(sigma)) unlist(sigma)
  if(is.list(Xo)) unlist(Xo)


	levels(regimes) <- 1:n_regimes
	regimes <- as.integer(regimes)-1  # convert to C-style indexing


	o<- .C("fit_model",
		as.double(Xo),
		as.double(alpha),
		as.double(theta),
		as.double(sigma),
		as.integer(regimes),
		as.integer(ancestor),
		as.double(branch_length),
		as.double(data),
		as.integer(n_nodes),
		as.integer(n_regimes),
		double(1),
		as.integer(use_siman) #use the simulated annealing approach 
	  )

	output <- list(data=dataIn, tree=tree, regimes=regimesIn, loglik=o[[11]],
              Xo = o[[1]], alpha = o[[2]], theta =  o[[3]], sigma = o[[4]])  
	class(output) <- "wrighttree"
	output
}






