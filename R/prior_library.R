

flat_prior_creator <- function(model_spec, n_regimes){
  indices <- get_indices(model_spec, n_regimes)
  prior<- function(pars){
    sum(2*pars[indices$alpha_i]/pars[indices$sigma_i]^2)
  }
}

general_prior <- function(pars){
# Prior for regime-specific alpha, theta, sigma
  n_regimes <- (length(pars) - 1)/3
## Is this really how you combine priors?  
  out <- dnorm(pars[1], 0, 1000, log=TRUE) +  # Xo prior -- 
  sum(dexp(pars[2:(1+n_regimes)], 1, log=TRUE)) +  # alpha prior -- 
  sum(dexp(pars[(2+n_regimes):(1+2*n_regimes)], 1, log=TRUE)) + # sigma prior -- 
  sum(dnorm(pars[(2+2*n_regimes):(1+3*n_regimes)], 0, 1000, log=TRUE)) # theta prior 
  out
}


# 
rc_prior <- function(pars){
# Prior for global theta, global sigma, regime-specific alpha
  n_regimes <- length(pars) - 3  

  out <- dnorm(pars[1], 0, 1000, log=TRUE) +  # Xo prior -- 
  sum(dexp(pars[2:(1+n_regimes)], 1, log=TRUE)) +  # alpha prior -- 
  dexp(pars[2+n_regimes], 1, log=TRUE) + # sigma prior -- 
  dnorm(pars[3+n_regimes], 0, 1000, log=TRUE) # theta prior 
  out
}
