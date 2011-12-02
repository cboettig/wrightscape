### Defines generics if they don't exist


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


#' An internal generic function
#' @keywords internal
getParameters <- function(x, y, ...) UseMethod("getParameters")

#' An internal S3 method to get the parameters from the wrightscape object
getParameters.wrighttree <- function(ws){
    c(alpha=ws$alpha, theta=ws$theta, sigma=ws$sigma, Xo=ws$Xo)
#, converge=ws$convergence) 
}

#' An internal S3 method to get the likelihood of a wrightscape object
loglik.wrighttree <- function(ws) ws$loglik



#' Simulate data under a multiOU model 
#' @param ws a multiOU model object (has elements tree, regimes, Xo, 
#'   alpha, theta, sigma).
#' @param ... additional parameters passed to internal simulate function
#' currently only "seed" is supported, gives random number generator seed
#' @return returns trait data (in ouch data format) from a simulate
#' @method simulate multiOU
#' @S3method simulate multiOU
#' @export
simulate.multiOU <- simulate.wrighttree


#' updates a multiOU model parameter estimation on new data
#' @param model a multiOU model
#' @param data trait data on which to base the update
#' @method update multiOU
#' @S3method update multiOU
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


#' Internal use mostly extending the method to multiOU class
#' exported for use by PMC
#' @param a multiOU class object (model fit)
#' @return extracts the log likelihood value
#' @method loglik multiOU
#' @S3method loglik multiOU
#' @export
loglik.multiOU <- loglik.wrighttree

#' Internal, extends method to mulitOU class
#' @param a multiOU class object
#' @return a vector of all the model parameter estimates
#' @method getParameters multiOU
#' @S3method getParameters multiOU
#' @export
getParameters.multiOU <- getParameters.wrighttree



