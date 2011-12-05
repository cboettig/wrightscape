\name{custom_multiType}
\alias{custom_multiType}
\title{get the likelihood of a custom model
custom_multiType allows the user to indicate which regimes will share
which parameters, making all submodels of the global model possible.}
\usage{
  custom_multiType(data, tree, regimes, par, indices,
    make_bm = NA, ...)
}
\arguments{
  \item{data}{the trait data}

  \item{tree}{the phylogenetic tree in ouch format}

  \item{regimes}{the regimes in ouch format}

  \item{par}{a list of initial parameters, structured
  according to indices (see details)}

  \item{indices}{a list with vectors alpha_i, sigma_i, and
  theta_i.  Each vector contains the index number of the
  parameter to be used for that regime -- repeated index
  numbers indicate a parameter that is the same for both
  regimes.  The entries in the vector are in order of the
  regimes and must be of length n_regimes.}

  \item{make_bm}{the regimes which should have their alpha
  value set to BM. the default value of NA will avoid
  forcing any to zero}

  \item{...}{extra options that are passed to the optimizer
  routine (optim).}
}
\description{
  get the likelihood of a custom model custom_multiType
  allows the user to indicate which regimes will share
  which parameters, making all submodels of the global
  model possible.
}
\details{
  Consider this example: Regimes are numbered 1 to 7.  1-6
  follow modelspec global alpha/sigma, indep theta and
  regime 7 follows modelspec alpha=0, theta, shares a
  global sigma.  Also imagine a parameter list, par =
  (alpha, sigma, theta1,2,3,4,5,6,7).  Then: indices =
  list(alpha_i = rep(1,7), sigma_i = rep(2,7), theta_i =
  3:(7+2)) make_bm <- c(7) # makes the 7th alpha value set
  to zero always
}
\examples{
## Define ouch's "hansen" function the custom way
## (equiv to using multiTypeOU withmodel_spec global gobal indep)
require(wrightscape)
require(ouch)
data(bimac)
tree <- with(bimac,ouchtree(node,ancestor,time/max(time),species))

par <- c(1, 1, 3, 3, 3) # init guesses for alpha, sigma, theta_1, theta_2, theta_3
# parameter mapping is created by indices:
indices <- list(alpha_i = c(1,1,1), sigma_i=c(2,2,2), theta_i = c(3,4,5))
ou3 <- custom_multiType(log(bimac[['size']]), tree, bimac[['OU.LP']], par, indices)
}

