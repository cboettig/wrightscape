\name{evaluate_likelihood}
\alias{evaluate_likelihood}
\title{get the likelihood of the specified model using the specified parameters}
\usage{
  evaluate_likelihood(data, tree, regimes,
    model_spec = list(alpha = "indep", sigma = "indep", theta = "indep"),
    Xo, alpha, sigma, theta)
}
\arguments{
  \item{data}{the trait data}

  \item{tree}{the phylogenetic tree in ouch format}

  \item{regimes}{the regimes in ouch format}

  \item{model_spec}{a list that specifies the model, see
  details}

  \item{Xo}{root state}

  \item{alpha}{a vector of length n_regimes if indep in
  model, or a scalar}

  \item{sigma}{a vector of length n_regimes if indep in
  model, or a scalar}

  \item{theta}{a vector of length n_regimes if indep in
  model, or a scalar}
}
\value{
  log likelihood
}
\description{
  get the likelihood of the specified model using the
  specified parameters
}
\details{
  the general model form is specified by model_spec list.
  This specifies which parameters out of alpha, theta, and
  sigma are independently estimated on each regime, kept
  global across regimes, or, in the case of alpha, fixed to
  zero (to give purely Brownian behavior).  i.e. ouch model
  is equivalent to: list(alpha="global", sigma="global",
  theta="indep"), while the brownie model is equivalent to
  list(alpha="fixed", sigma="indep", theta="global")
}

