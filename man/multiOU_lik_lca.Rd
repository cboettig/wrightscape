\name{multiOU_lik_lca}
\alias{multiOU_lik_lca}
\title{compute the likelihood by passing data to the C level function, wrightscape}
\usage{
  multiOU_lik_lca(data, tree, regimes, alpha = NULL,
    sigma = NULL, theta = NULL, Xo = NULL, lca)
}
\arguments{
  \item{data}{- ouch-style data}

  \item{tree}{- ouch-tree}

  \item{regimes}{- painting of selective regimes, as in
  ouch}

  \item{alpha}{- (vector length n_regimes) gives strength
  of selection in each regime}

  \item{sigma}{- (vector length n_regimes) gives
  diversification rate in each regime}

  \item{theta}{- (vector length n_regimes) gives optimum
  trait in each regime}

  \item{Xo}{- root value}

  \item{lca}{- least common ancestor matrix, from lca_calc
  fun}
}
\value{
  the log likelihood at the given parameter values
}
\description{
  compute the likelihood by passing data to the C level
  function, wrightscape
}
\details{
  A version of the multitype OU likelihood function that
  accepts the least common ancestor matrix as a parameter.
  This may increase computational speed, since the
  calculation needs to be done only once and is rather slow
  as implemented .C calls should do some error checking on
  the length of inputs maybe, to avoid crashes when given
  inappropriate calls
}

