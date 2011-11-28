\name{simulate.multiOU}
\alias{simulate.multiOU}
\title{Simulate data under a multiOU model}
\usage{
  \method{simulate}{multiOU} (ws, ...)
}
\arguments{
  \item{ws}{a multiOU model object (has elements tree,
  regimes, Xo, alpha, theta, sigma).}

  \item{...}{additional parameters passed to internal
  simulate function currently only "seed" is supported,
  gives random number generator seed}
}
\value{
  returns trait data (in ouch data format) from a simulate
}
\description{
  Simulate data under a multiOU model
}

