\name{summary.multiOU}
\alias{summary.multiOU}
\title{Summarize a multiOU fit's outputs, including bootstraps, if provided}
\usage{
  \method{summary}{multiOU} (modelfit, bootstrap = NULL,
    silent = FALSE)
}
\arguments{
  \item{modelfit}{a multiOU class model fit (from
  multiTypeOU)}

  \item{bootstrap}{the bootstraps of the object. Optional,
  otherwise just gives fit.}

  \item{silent}{a logical indicating if summary info should
  be printed to terminal}
}
\value{
  a list containing parameters estimated, the summary of
  the bootstrap, and the bootstrap object, if provided.
  Otherwise just parameters and a message that bootstrap
  wasn't provided.
}
\description{
  Summarize a multiOU fit's outputs, including bootstraps,
  if provided
}
\examples{
data(parrotfish)
 alphas <- multiTypeOU(data=dat, tree=tree, regimes = intramandibular,
  model_spec = list(alpha = "indep",sigma = "global", theta = "indep"))
 boots <- replicate(3, bootstrap(alphas))
 summary(alphas, boots)
}

