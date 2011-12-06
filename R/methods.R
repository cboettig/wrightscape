### Defines generics if they don't exist

#' simulate the multiOU model
#' @return returns trait data (in ouch data format) from a simulate
#' 
# @rdname simulate
# @method simulate multiOU
# @S3method simulate multiOU
#' @export
simulate.multiOU <-  function(ws, ...){
# @param ws a multiOU model object (has elements tree, regimes, Xo, 
#   alpha, theta, sigma).
# @param ... additional parameters passed to internal simulate function
# currently only "seed" is supported, gives random number generator seed
	output <- simulate_wrightscape(tree=ws$tree, regimes=ws$regimes,
                                       Xo=ws$Xo, alpha=ws$alpha,
                                       theta=ws$theta, sigma=ws$sigma, ...)
  output$rep.1
}

#' An internal generic function
#' @title get parameters of a fitted object
#' @param x a model fit 
#' @param ... additional arguments
#' @export
getParameters <- function(x, ...) UseMethod("getParameters")

#' get parameters for the multiOU model
#' @return a vector of parameters
#' 
# @rdname getParameters
# @method getParameters multiOU
# @S3method getParameters multiOU
#' @export
getParameters.multiOU <- function(ws){
    c(alpha=ws$alpha, theta=ws$theta, sigma=ws$sigma, Xo=ws$Xo, converge=ws$convergence) 
}

#' get the log likelihood
#' @return log likelihood
#'
# @rdname loglik
# @method loglik multiOU
# @S3method loglik multiOU
#' @export
loglik.multiOU <- function(ws) ws$loglik

#' update the model fit 
#' @return updated model fit 
#'
# @rdname update
# @method update multiOU
# @S3method update multiOU
#' @export
update.multiOU <- function(model, data){
# Define the update method for this style of fit.  
# Simulate, GetParNames, loglik methods all defined in wrightscape.R
  do.call(multiTypeOU, c(list(data=data, tree=model$tree,
                                   regimes=model$regimes,
                                   model_spec=model$model_spec,
                                   Xo=model$Xo, alpha=model$alpha, 
                                  sigma=model$sigma, theta=model$theta),
                            model$opts))
}





## old class, wrighttree, deprecated 

# An internnal S3 method, updates with new data
update.wrighttree <- function(ws, data){
	wrightscape(data=data, tree=ws$tree, regimes=ws$regimes, alpha=ws$alpha,
              	    sigma=ws$sigma, theta=ws$theta, Xo=ws$Xo)
}

#' An internal S3 method, simulates the data of a wrightscape wrighttree object
simulate.wrighttree <- function(ws, ...){
	output <- simulate_wrightscape(tree=ws$tree, regimes=ws$regimes,
                                       Xo=ws$Xo, alpha=ws$alpha,
                                       theta=ws$theta, sigma=ws$sigma, ...)
  output$rep.1
}

#' An internal S3 method to get the likelihood of a wrightscape object
loglik.wrighttree <- function(ws) ws$loglik



#' An internal S3 method to get the parameters from the wrightscape object
getParameters.wrighttree <- function(ws){
    c(alpha=ws$alpha, theta=ws$theta, sigma=ws$sigma, Xo=ws$Xo, converge=ws$convergence) 
}


