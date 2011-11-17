\name{format_data}
\alias{format_data}
\title{format data in ape format into ouch format}
\usage{
  format_data(tree, traits, species_names = NULL,
    regimes = NULL)
}
\arguments{
  \item{tree}{a phylogenetic tree of class "phylo", ape
  format}

  \item{traits}{a numeric with trait values, or a matrix or
  data frame of traits, rownames matching species or handed
  in}

  \item{species_names}{in the order of entries in traits,
  if not given in rownames.}

  \item{regimes}{the column in traits containing regime
  labels}
}
\value{
  the ouch-formatted tree, traits, and regimes
}
\description{
  format data in ape format into ouch format
}
\details{
  Should become an internal function to handle data
  conversion to ape-type
}
\seealso{
  \code{\link{convert}} to toggle between formats,
  including regime paintings
}

