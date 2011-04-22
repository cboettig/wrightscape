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

    ## calculate the likelihood
	o<- .C("calc_lik", as.double(Xo), as.double(alpha), as.double(theta),
            as.double(sigma), as.integer(regimes), as.integer(ancestor),
            as.double(branch_length), as.double(data), as.integer(n_nodes),
            as.integer(lca[[4]]))

	output <- list(data=dataIn, tree=tree, regimes=regimesIn, loglik=o[[11]], Xo = o[[1]], alpha = o[[2]], theta =  o[[3]], sigma = o[[4]]  )  
	class(output) <- "wrighttree"
	output
}
