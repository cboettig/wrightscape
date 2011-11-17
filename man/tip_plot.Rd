\name{tip_plot}
\alias{tip_plot}
\title{creates a color-coded plot of how well each tip fits it's predicted val}
\usage{
  tip_plot(modelfit, n = 100, ...)
}
\arguments{
  \item{modelfit}{a fit from multiTypeOU (class multiOU)}

  \item{n}{number of replicates to use in simulation}

  \item{...}{additional arguments to plot.phylo}
}
\description{
  creates a color-coded plot of how well each tip fits it's
  predicted val
}
\details{
  this function is a proof of principle that could be made
  more generic Note that this function also illustrates
  tree conversion with regime data.
}

