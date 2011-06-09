# Define generics if they dont exist
getParameters <- function(x, y, ...) UseMethod("getParameters")


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

# update with new data
update.wrighttree <- function(ws, data){
	wrightscape(data=data, tree=ws$tree, regimes=ws$regimes, alpha=ws$alpha,
              sigma=ws$sigma, theta=ws$theta, Xo=ws$Xo)
}


simulate.wrighttree <- function(ws, ...){
	output <- simulate_wrightscape(tree=ws$tree, regimes=ws$regimes,
                                   Xo=ws$Xo, alpha=ws$alpha,
                                   theta=ws$theta, sigma=ws$sigma, ...)
  output$rep.1
}
loglik.wrighttree <- function(ws) ws$loglik

getParameters.wrighttree <- function(ws){
    c(alpha=ws$alpha, theta=ws$theta, sigma=ws$sigma, Xo=ws$Xo)
#, converge=ws$convergence) 
}

simulate.multiOU <- simulate.wrighttree
loglik.multiOU <- loglik.wrighttree
getParameters.multiOU <- getParameters.wrighttree



