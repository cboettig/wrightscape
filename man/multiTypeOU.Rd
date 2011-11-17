\name{multiTypeOU}
\alias{multiTypeOU}
\title{Fits the generic multitype OU.}
\usage{
  multiTypeOU(data, tree, regimes,
    model_spec = list(alpha = "indep", sigma = "indep", theta = "indep"),
    Xo = NULL, alpha = 1, sigma = 1, theta = NULL, ...)
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
  Fits the generic multitype OU.
}
\details{
  Submodels such as brownie, and other unique models, can
  be created by specifiying how parameters are treated
  under model_spec get the likelihood of the specified
  model using the specified parameters
}

