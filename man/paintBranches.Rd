\name{paintBranches}
\alias{paintBranches}
\title{Paints branches with regimes changing at nodes specified}
\usage{
  paintBranches(regimeShiftNodes, tree,
    regimeTitles = NULL)
}
\arguments{
  \item{tree}{OUCH-style (S4) tree}

  \item{regimeShiftNodes}{a vector of nodes or a list of
  taxa defining nodes at which selective regimes shift:
  root may be included--if not present, will be added--but
  tips are meaningless in this context}

  \item{regimeTitles}{a vector of titles for the regimes
  that begin at the root and at the nodes indicated in
  "regimeShiftNodes", in order of description in
  "regimeShiftNodes", except that the root is listed first
  in "regimeTitles" but not at all in "regimeShiftNodes"...
  defaults to "node[x]regime}
}
\value{
  a vector of regimes that can be handed to hansen
}
\description{
  Paints branches with regimes changing at nodes specified
}

