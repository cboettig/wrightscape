\name{custom.llik.closure}
\alias{custom.llik.closure}
\title{returns a custom-built likelihood function}
\usage{
  custom.llik.closure(data, tree, regimes, indices,
    make_bm)
}
\arguments{
  \item{data}{the trait data}

  \item{tree}{the phylogenetic tree in ouch format}

  \item{regimes}{the regimes in ouch format}

  \item{indices}{a list with vectors alpha_i, sigma_i, and
  theta_i.  Each vector contains the index number of the
  parameter to be used for that regime -- repeated index
  numbers indicate a parameter that is the same for both
  regimes.  The entries in the vector are in order of the
  regimes and must be of length n_regimes.}

  \item{make_bm}{the regimes which should have their alpha
  value set to BM.}
}
\value{
  the negative log likelihood of the parameters given the
  data.
}
\description{
  returns a custom-built likelihood function
}
\details{
  Consider this example: Regimes are numbered 1 to 7.  1-6
  follow modelspec global alpha/sigma, indep theta and
  regime 7 follows modelspec alpha=0, theta, shares global
  sigma imagine a parameter list: (alpha, sigma,
  theta1,2,3,4,5,6,7) indices = list(alpha_i = rep(1,7),
  sigma_i = rep(2,7), theta_i = 3:(7+2)) make_bm <- c(7) #
  makes the 7th alpha value set to zero always
}

