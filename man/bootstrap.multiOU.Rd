\name{bootstrap.multiOU}
\alias{bootstrap.multiOU}
\title{function to bootstrap a multiOU object}
\usage{
  \method{bootstrap}{multiOU} (modelfit)
}
\arguments{
  \item{modelfit}{a multiOU object, (output of
  multiTypeOU() fn)}
}
\value{
  a matrix of bootstrap values. (replicates fixed values)
}
\description{
  function to bootstrap a multiOU object
}
\examples{
data(parrotfish)
 alphas <- multiTypeOU(data=dat, tree=tree, regimes = intramandibular,
  model_spec = list(alpha = "indep",sigma = "global", theta = "indep"))
 boots <- replicate(3, bootstrap(alphas))
 summary(alphas, boots)
}

