#multiou can try and take lca as a parameter option rather than calculating each time, for efficiency

wright <- function(data, tree, regimes, Xo=NULL, alpha=1, sigma=1){
    # intialize a parameter vector to optimize: 
    # Xo, alpha, sigma, and the n_regime thetas
    n_regimes <- length(levels(regimes))
    par <- numeric(1+3*n_regimes)

    # Some starting conditions
    if(is.null(Xo)) Xo <- mean(data, na.rm=TRUE) 
    par[1] <- Xo
    par[2:(1+n_regimes)] <- rep(alpha, n_regimes)
    par[(2+n_regimes):(1+2*n_regimes)] <- rep(sigma, n_regimes)
    par[(2+2*n_regimes):(1+3*n_regimes)] <- rep(Xo, n_regimes)

    # Likelihood as a function of optimizable parameters
    f <- function(par){
        Xo <- par[1]
        alpha <- par[2:(1+n_regimes)]
        sigma <- par[(2+n_regimes):(1+2*n_regimes)]
        theta <- par[(2+2*n_regimes):(1+3*n_regimes)] 
        o<-multiOU_lik(data, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo)
        -o$loglik
    }
    optim(par,f) 
}




# OUCH
ouch <- function(data, tree, regimes, Xo=NULL, alpha=1, sigma=1){
    # intialize a parameter vector to optimize: 
    # Xo, alpha, sigma, and the n_regime thetas
    n_regimes <- length(levels(regimes))
    par <- numeric(3+n_regimes)

    # Some starting conditions
    if(is.null(Xo)) Xo <- mean(data, na.rm=TRUE) 
    par[1] <- Xo
    par[2] <- alpha
    par[3] <- sigma
    par[-c(1:3)] <- rep(Xo, n_regimes)

    # Likelihood as a function of optimizable parameters
    f <- function(par){
        Xo <- par[1]
        alpha <- rep(par[2], n_regimes)
        theta <- par[-c(1:3)]
        sigma <- rep(par[3], n_regimes) # everything else
        o <- multiOU_lik(data, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo)
        -o$loglik
    }
    optim(par,f) 
}

# Brownie
brownie <- function(data, tree, regimes, sigma=1){ 

    # intialize a parameter vector to optimize: 
    # Xo, followed by the n_regime sigmas
    n_regimes <- length(levels(regimes))
    par <- numeric(1+n_regimes)

    # Some starting conditions
    par[1] <- mean(data, na.rm=TRUE) #Xo
    par[-1] <- rep(sigma, n_regimes) # sigmas

    # Likelihood as a function of optimizable parameters
    f <- function(par){
        Xo <- par[1]
        sigma <- par[-1] # everything else
        alpha <- rep(1e-12, n_regimes) 
        theta <- rep(Xo, n_regimes)
        o <- multiOU_lik(data, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo)
        -o$loglik
    }
    optim(par,f) 
}

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

	output <- list(data=dataIn, tree=tree, regimes=regimesIn, loglik=o[[11]], Xo = o[[1]], alpha = o[[2]], theta =  o[[3]], sigma = o[[4]]  )  
	class(output) <- "wrighttree"
	output
}
