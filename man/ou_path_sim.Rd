\name{ou_path_sim}
\alias{ou_path_sim}
\title{OU path simulator}
\usage{
  ou_path_sim(x0 = 0, T = 1, alpha = 4, sigma = 1,
    theta = 0, pts = 100, reps = 40)
}
\arguments{
  \item{x0}{starting position}

  \item{T}{total time}

  \item{alpha}{selection/constraint strength}

  \item{sigma}{noise parameter}

  \item{theta}{optimum value}

  \item{pts}{number of points per rep}

  \item{reps}{number of replicates}
}
\value{
  a melted data object ready for ggploting
}
\description{
  OU path simulator
}

