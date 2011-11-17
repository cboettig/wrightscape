# File: bootstrap.R
# Author: Carl Boettiger <cboettig@gmail.com>
# Date: 2011-11-16
# License: BSD

#' @title bootstrap: a function to compute bootstrap replicates
#' define the generic S3 bootstrap method
bootstrap <- function(x, ...) UseMethod("bootstrap")

#' function to bootstrap a multiOU object
#' @param modelfit a multiOU object, (output of multiTypeOU() fn)
#' @returns a matrix of bootstrap values. (replicates fixed values)
#' @examples
#'  data(parrotfish)
#'  alphas <- multiTypeOU(data=dat, tree=tree, regimes = intramandibular, 
#'   model_spec = list(alpha = "indep",sigma = "global", theta = "indep"))
#'  boots <- replicate(3, bootstrap(alphas))
#'  summary(alphas, boots)
#' @S3method bootstrap multiOU
#' @method bootstrap multiOU
#' @export
bootstrap.multiOU <- function(modelfit){
      dat <- simulate(modelfit) 
      out <- update(modelfit, dat)

      n <- length(levels(modelfit$regimes))
      Xo <- rep(out$Xo,n) 
      loglik <- rep(out$loglik, n)
      pars <- rep(NA, 5 * n)
      if(out$convergence == 0) # only return values if successful
        pars <- c(out$alpha, out$sigma, out$theta, Xo, loglik)
      pars
    }

#' Summarize a multiOU fit's outputs, including bootstraps, if provided
#' @param modelfit a multiOU class model fit (from multiTypeOU)
#' @param bootstrap the bootstraps of the object. Optional, otherwise just gives fit.
#' @param silent a logical indicating if summary info should be printed to terminal
#' @return a list containing parameters estimated, the summary of the bootstrap, and
#'  the bootstrap object, if provided. Otherwise just parameters and a message that
#'  bootstrap wasn't provided. 
#' @examples
#'  data(parrotfish)
#'  alphas <- multiTypeOU(data=dat, tree=tree, regimes = intramandibular, 
#'   model_spec = list(alpha = "indep",sigma = "global", theta = "indep"))
#'  boots <- replicate(3, bootstrap(alphas))
#'  summary(alphas, boots)
#' @S3method 
summary.multiOU <- function(modelfit, bootstrap = NULL, silent = FALSE){
  est <- rbind(alpha = modelfit$alpha, sigma = modelfit$sigma,
               theta = modelfit$theta)
  if(is.null(bootstrap)){
    SE <- "Bootstraps not provided, run bootstrap(modelfit) to get SE"
    bootstrap <- "Not provided"
  } else if(is.matrix(bootstrap)){

    SE <- sapply(1:dim(bootstrap)[1], function(i) sd(bootstrap[i,], na.rm=T) )
    SE <- t(matrix(SE, nrow = length(levels(modelfit$regimes))))
    rownames(SE) = c("alpha", "sigma", "theta", "Xo", "loglik") 
    colnames(SE) = levels(modelfit$regimes)
    colnames(est) = levels(modelfit$regimes)
  } else {
    warning("bootstrap type not recognized")
  }
  if(!silent){
    print("Parameter Estimates")
    print(est)
    print("Standard Error Estimates")
    print(SE)
    print(paste("log likelihood:", modelfit$loglik, "Root state Xo:", modelfit$Xo))
  }
  list(Param.est = est, Param.SE = SE, bootstraps=bootstrap)
}



