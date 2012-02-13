\name{release_path_sim}
\alias{release_path_sim}
\title{Release of constraint path simulator}
\usage{
  release_path_sim(x0 = 0, T = 1, alpha = 4, sigma = 1,
    theta = 0, pts = 100, reps = 40, release_frac = 0.5)
}
\arguments{
  \item{x0}{starting position}

  \item{T}{total time}

  \item{alpha}{selection/constraint strength}

  \item{sigma}{noise parameter}

  \item{theta}{optimum value}

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

