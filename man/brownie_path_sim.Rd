\name{brownie_path_sim}
\alias{brownie_path_sim}
\title{Release of constraint path simulator}
\usage{
  brownie_path_sim(x0 = 0, T = 1, sigma = 1/sqrt(8),
    sigma2 = 1, pts = 100, reps = 40, release_frac = 0.5)
}
\arguments{
  \item{x0}{starting position}

  \item{T}{total time}

  \item{sigma}{BM rate at start}

  \item{sigma2}{BM rate after shift}

  \item{pts}{number of points per rep}

  \item{reps}{number of replicates}

  \item{release_frac}{fraction of time after which shift
  occurs}
}
\value{
  a melted data object ready for ggploting
}
\description{
  Release of constraint path simulator
}

