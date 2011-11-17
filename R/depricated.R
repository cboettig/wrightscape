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











## These should be part of an independent phylogenetic bootstrapping library (or at least file)
## Consider extending to some other functions, such as ape's ace fn, geiger's ancestral states, etc

LR_bootstrap <- function(true_model, test_model, nboot = 200){
# Bootstraps the likelihood ratio statistic using boot function
# Args:
#		true_model -- is used to generated the simulated data.  Must be a fitted hansentree or browntree
#		test_model -- is another fitted model whose likeihood will also be evaluated on the data
#		nboot -- is the number of bootstrap replicates to do.  Defaults to 200
#	Returns:
#		boot object that can be fed into boot.ci to 
#			generate confidence intervals, etc using the boot package

	get_loglik <- function(model){
		if(is(model, "ouchtree") ) loglik = model@loglik 
		else loglik = model$loglik
		loglik
	}
	get_data <- function(model){
		if(is(model, "ouchtree") ) data = model@data 
		else data = model$data
		data
	}

	orig_diff <- -2*( get_loglik(true_model) - get_loglik(test_model))
	orig_data <- get_data(true_model)

	statisticfn <- function(data, ...){
	# function required by boot fn that will be boostrapped
	# Args: 
	#		data is data generated by a simulation, 
	#		... is the test_model
	# Returns 
	#		the likelihood ratio statistic: 2*(log(null) - log(test) )
		test <- update(test_model, data=data)
		true <- update(true_model,data=data)
		-2*(get_loglik(true) - get_loglik(test))
	}

	rangendat <- function(d, p){
	# simulate using model specified as mle, the true model
		out <- simulate(p)
		out$rep.1
	}

	boot.out <- boot(	data=orig_data, 
						statistic=statisticfn, 
						R=nboot, 
						sim="parametric", 
						ran.gen=rangendat, 
						mle=true_model, 
						object=test_model)
}


fast_boot <- function(model, nboot=200, cpus=1){

	require(snowfall)
	sfInit(parallel=TRUE, cpus=cpus)
	sfExportAll()
	sfLibrary(wrightscape)

	fits <- sfLapply(1:nboot, function(i) update(model, data=simulate(model)$rep.1 ) )
	if(is(fits[[1]], "wrighttree") ){
		n_regimes <- length( fits[[1]]$sigma )
		n_pars <- 3*n_regimes+2
		regime_names <- levels(fits[[1]]$regimes[[1]])
		alpha_names <- sfSapply(1:n_regimes, function(i) paste("alpha.", regime_names[i]) )
		sigma_names <- sfSapply(1:n_regimes, function(i) paste("sigma.", regime_names[i]) )
		theta_names <- sfSapply(1:n_regimes, function(i) paste("theta.", regime_names[i]) )

		X <- sapply(1:nboot, function(i)  c(fits[[i]]$loglik, fits[[i]]$Xo, fits[[i]]$alpha, fits[[i]]$sigma, fits[[i]]$theta ) )
		rownames(X) <- c("loglik", "Xo", alpha_names, sigma_names, theta_names) 
	}
	out <- list(bootstrap_values = X, model=model, nboot=nboot)
	class(out) <- "wrightboot"
	out

}


plot.wrightboot <- function(input, CHECK_OUTLIERS=FALSE){
	object <- input$bootstrap_values
	par(mfrow=c(1,3) )
	n_regimes <- (dim(object)[1]-2)/3
	alphas <- 3:(2+n_regimes)
	sigmas <- (3+n_regimes):(2*n_regimes+2)
	thetas <- (3+2*n_regimes):(3*n_regimes+2)
	nboot <- dim(object)[2]
	outliers <- numeric(nboot)

	xlim <- c(0, 3*median(object[alphas,]) )
	ylim <- c(0, max(sapply(alphas, function(i) max(density(object[i,])$y))))

	plot(density(object[alphas[1], ]), xlim=xlim, ylim=ylim, xlab="Alpha values", type='n', main="", cex.lab=1.6, cex.axis = 1.6)
	k <- 1
	for(i in alphas){
		if(CHECK_OUTLIERS) outliers <- object[i,] > xlim[2]
		if( sum(outliers) > 0 ) print(paste(sum(outliers), " outliers in alpha ", k))
		lines(density(object[i,!outliers]), lwd = 3, lty=k)
		k <- k+1
	}

	xlim <- c(0, 3*median(object[sigmas,] ) )
	ylim <- c(0, max(sapply(sigmas, function(i) max(density(object[i,])$y))))
	plot(density(object[sigmas[1], ]), xlim=xlim, ylim=ylim, xlab="Sigma values", type='n', main="", cex.lab=1.6, cex.axis = 1.6)
	k <- 1
	for(i in sigmas){
		if(CHECK_OUTLIERS) outliers <- object[i,] > xlim[2]
		if( sum(outliers) > 0 ) print(paste(sum(outliers), " outliers in sigma ", k))
		lines(density(object[i,!outliers]), lwd = 3, lty=k)
		k <- k+1
	}

	xlim <- c(0, 3*median(object[thetas,] ) )
	ylim <- c(0, max(sapply(thetas, function(i) max(density(object[i,])$y))))
	plot(density(object[thetas[1], ]), xlim=xlim, ylim=ylim, xlab="Theta values", type='n', main="", cex.lab=1.6, cex.axis = 1.6)
	k <- 1
	for(i in thetas){
		if(CHECK_OUTLIERS) outliers <- object[i,] > xlim[2]
		if( sum(outliers) > 0 ) print(paste(sum(outliers), " outliers in theta ", k))
		lines(density(object[i,!outliers]), lwd = 3, lty=k)
		k <- k+1
	}
	legend("topright", levels(input$model$regimes), lty=1:length(levels(input$model$regimes)), lwd=2) 

}


bootstrap.wrighttree <- function(model, nboot = 200, fit=TRUE)
{
# Bootstraps the likelihood ratio statistic using boot function.  Should give bootstraps for all parameters!!!!
# Args:
#		model -- is used to generated the simulated data.  Must be a fitted hansentree or browntree
#		nboot -- is the number of bootstrap replicates to do.  Defaults to 200
#	Returns:
#		boot object -- can be fed into boot.ci to generate confidence intervals, etc 
#			using the boot package

	simdata <- simulate(model)
	refit_model <- update(model, data=simdata)
}

choose_model <- function(model_list, nboot=200, cpus=1){
	require(snowfall)
	sfInit(parallel=TRUE, cpus=cpus)
	sfExportAll()
	sfLibrary(wrightscape)

	LR <- sfLapply( 1:(length(model_list)-1),
					 function(i) LR_bootstrap( model_list[[i]], model_list[[i+1]], nboot )
		  		   )
	p_vals <- sfSapply( 1:(length(model_list)-1),
			function(i)  sum( LR[[i]]$t < LR[[i]]$t0 )/length(LR[[i]]$t)
		  )
	print(p_vals)
	LR	
}


pretty_plot <- function(LR, main=""){
	xlim = 1.1*c(min( LR$t, LR$t0), max( LR$t, LR$t0) )
	hist(LR$t, col="lightblue", border="white", xlab="Likelihood Ratio", cex.axis=1.6, cex.lab=1.6, main=main, xlim=xlim)
	abline(v=LR$t0, lwd=4, lty=2, col="darkblue")
	p_val <- 1-sum(LR$t < LR$t0)/length(LR$t)  
	text(LR$t0, 0.5*par()$yaxp[2], paste("p = ", round(p_val, digits=3)), cex=1.6)   # halfway up the vert line
}

# plot the wrightscape tree using the ouch plotting function
plot.wrighttree <- function(object)
{
	plot(object$tree, regimes=object$regimes)
}


