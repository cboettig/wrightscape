\name{paint_phy}
\alias{paint_phy}
\title{returns a phylogenetic tree painted for ouwie and data}
\usage{
  paint_phy(phy, data, clades, show_plot = TRUE)
}
\arguments{
  \item{phy}{a phylogenetic tree of class phylo (ape)}

  \item{clades:}{a list of clades defined by species pairs,
  e.g. list(c(sp1, sp2), c(sp3, sp4)), or just a vector
  c(sp1, sp2) the descendents of the MRCA of sp1 and sp2
  define the regime later clades are painted over earlier
  ones, so should be in order of nesting or otherwise
  non-overlapping.}

  \item{data}{the data.frame of characters (see data in
  ?treedata)}

  \item{show_plot}{logical, should I plot the resulting
  tree}
}
\value{
  a list with the phy, updated to have regimes at node
  labels (for OUwie) and the data matrix formatted with
  regimes, for OUwie, and a list of colors corresponding to
  edges, that can be used for plotting.
}
\description{
  returns a phylogenetic tree painted for ouwie and data
}
\examples{
data(labrids)
myphy <- paint_phy(ape.phy, ape.dat,
                   list(c("Bolbometopon_muricatum", "Sparisoma_radians"),
                   c("Chlorurus_sordidus", "Hipposcarus_longiceps")))
plot(myphy$phy, edge.color=myphy$colors, type="fan",
show.tip.label=FALSE, edge.width=2)
}

