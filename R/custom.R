#' returns a custom-built likelihood function
#' @param data the trait data
#' @param tree the phylogenetic tree in ouch format
#' @param regimes the regimes in ouch format
#' @param indices a list with vectors alpha_i, sigma_i, and theta_i.  Each vector
#' contains the index number of the parameter to be used for that regime -- repeated
#' index numbers indicate a parameter that is the same for both regimes.  The entries
#' in the vector are in order of the regimes and must be of length n_regimes. 
#' @param make_bm the regimes which should have their alpha value set to BM.  
#' @return the negative log likelihood of the parameters given the data. 
#' @details Consider this example:
#' Regimes are numbered 1 to 7.  1-6 follow modelspec global alpha/sigma, indep theta
#' and regime 7 follows modelspec alpha=0, theta, shares global sigma
#' imagine a parameter list: (alpha, sigma, theta1,2,3,4,5,6,7)
#' indices = list(alpha_i = rep(1,7), sigma_i = rep(2,7), theta_i = 3:(7+2)) 
#' make_bm <- c(7)  # makes the 7th alpha value set to zero always
custom.llik.closure <- function(data, tree, regimes, indices, make_bm){
  n_regimes <- length(levels(regimes))
  lca <- lca_calc(tree)
  f <- function(par){

    Xo <- par[indices$theta_i][match(regimes[1], levels(regimes))]  # assumes root theta
    alpha <- par[indices$alpha_i]
    alpha[make_bm] <- 1e-12
    sigma <- par[indices$sigma_i]
    theta <- par[indices$theta_i] 

    # avoid nonsense parameter estimates
    if (any(alpha < 0))
        llik <- -Inf
    else if (any(sigma<0))
        llik <- -Inf

    ## All set, here we go a-calculating:
    else 
        llik<-multiOU_lik_lca(data, tree, regimes, alpha=alpha,
                              sigma=sigma, theta=theta, Xo=Xo, lca)

    # more avoiding nonsense estimates
    if(is.nan(llik))
      llik <- -Inf
    # returns the negative log likelihood value 
    -llik
  }
  # returns a likelihood function of the parameters for the optimizer routine
  f
}


#' get the likelihood of a custom model
#' custom_multiType allows the user to indicate which regimes will share 
#' which parameters, making all submodels of the global model possible.  
#' @param data the trait data
#' @param tree the phylogenetic tree in ouch format
#' @param regimes the regimes in ouch format
#' @param par a list of initial parameters, structured according 
#' to indices (see details)
#' @param indices a list with vectors alpha_i, sigma_i, and theta_i.  Each vector
#' contains the index number of the parameter to be used for that regime -- repeated
#' index numbers indicate a parameter that is the same for both regimes.  The entries
#' in the vector are in order of the regimes and must be of length n_regimes. 
#' @param make_bm the regimes which should have their alpha value set to BM.
#' the default value of NA will avoid forcing any to zero 
#' @param ... extra options that are passed to the optimizer routine (optim).   
#' @details  Consider this example:
#' Regimes are numbered 1 to 7.  1-6 follow modelspec global alpha/sigma, 
#' indep theta and regime 7 follows modelspec alpha=0, theta, shares a 
#' global sigma.  Also imagine a parameter list,
#' par <- (alpha, sigma, theta1,2,3,4,5,6,7).  Then:
#' indices <- list(alpha_i = rep(1,7), sigma_i = rep(2,7), theta_i = 3:(7+2)) 
#' make_bm <- c(1)  # makes the 7th alpha value set to zero always
#' @useDynLib wrightscape
#' @import ouch
#' @import geiger
#' @examples 
#' ## Define ouch's "hansen" function the custom way 
#' ## (equiv to using multiTypeOU withmodel_spec global gobal indep)
#' require(wrightscape)
#' require(ouch)
#' data(bimac)
#' tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))
#' 
#' par <- c(1, 1, 3, 3, 3) # init guesses for alpha, sigma, theta_1, theta_2, theta_3 
#' # parameter mapping is created by indices:
#' indices <- list(alpha_i = c(1,1,1), sigma_i=c(2,2,2), theta_i = c(3,4,5))
#' ou3 <- custom_multiType(log(bimac[['size']]), tree, bimac[['OU.LP']], par, indices)
#' @export
custom_multiType <- function(data, tree, regimes, par, indices, make_bm=NA,   
                     ...){

  myCall <- match.call() # can be used to call the function as it was called

  # Create the submodel and do the optimization (all the action is here)
  f <- custom.llik.closure(data, tree, regimes, indices, make_bm)
  optim_output <- optim(par, f, ...)

  # extract parameter list
  out.par <- optim_output$par
  Xo <- out.par[indices$theta_i][match(regimes[1], levels(regimes))]  
  alpha <- out.par[indices$alpha_i]
  alpha[make_bm] <- 1e-12
  sigma <- out.par[indices$sigma_i]
  theta <- out.par[indices$theta_i] 

  names(alpha) <- levels(regimes)
  names(theta) <- levels(regimes)
  names(sigma) <- levels(regimes)

  output <- list(data=data, tree=tree, regimes=regimes, 
                 loglik=-optim_output$value, alpha=alpha, sigma=sigma, 
                 theta=theta, Xo = Xo, optim_output=optim_output, 
                 convergence=optim_output$convergence, 
                 indices=indices, make_bm=make_bm,
                 myCall = myCall)
  class(output) = "multiOU"
  output
}



