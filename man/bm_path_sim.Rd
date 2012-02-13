\name{bm_path_sim}
\alias{bm_path_sim}
\title{brownian motion path simulator}
\usage{
  bm_path_sim(x0 = 0, T = 1, sigma = 1, pts = 100,
    reps = 40)
}
\arguments{
  \item{x0}{starting position}

  \item{T}{total time}

  \item{sigma}{noise parameter}

  \item{pts}{number of points per rep}

  \item{reps}{number of replicates}
}
\value{
  a melted data object ready for ggploting
}
\description{
  brownian motion path simulator
}

