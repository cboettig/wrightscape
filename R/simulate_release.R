#' brownian motion path simulator
#'
#' @param x0 starting position
#' @param T total time
#' @param sigma noise parameter 
#' @param pts number of points per rep
#' @param reps number of replicates
#' @return a melted data object ready for ggploting 
#' @export
bm_path_sim <- function(x0 = 0, T = 1, sigma = 1, pts = 100, reps = 40){
  deltaT <- T/pts
  X <- replicate(reps,  x0 + cumsum(rnorm(pts, 0, sigma*sqrt(deltaT))))
  out <- melt(X)
  names(out) <- c("time", "rep", "value")
  out
}


#' OU path simulator
#'
#' @param x0 starting position
#' @param T total time
#' @param alpha selection/constraint strength 
#' @param sigma noise parameter 
#' @param theta optimum value
#' @param pts number of points per rep
#' @param reps number of replicates
#' @return a melted data object ready for ggploting 
#' @export
ou_path_sim <- function(x0 = 0, T = 1, alpha = 4, sigma = 1, theta = 0,
                        pts = 100, reps = 40){
  deltaT <- T/pts
  ou <- function(){
    x <- numeric(pts)
    x[1] <- x0
    for(i in 1:pts){
      mu <- x[i] * exp(-alpha * deltaT) + theta * (1 - exp(-alpha * deltaT))
      var <- sigma * sigma * (1 - exp(-2 * alpha * deltaT)) / (2 * alpha) 
      x[i + 1] <- rnorm(1, mu, sqrt(var))
    }
    x
  }
  X <- replicate(reps, ou())
  out <- melt(X)
  names(out) <- c("time", "rep", "value")
  out
}


#' Release of constraint path simulator
#'
#' @param x0 starting position
#' @param T total time
#' @param alpha selection/constraint strength 
#' @param sigma noise parameter 
#' @param theta optimum value
#' @param pts number of points per rep
#' @param reps number of replicates
#' @param release_frac fraction of time after which shift occurs
#' @return a melted data object ready for ggploting 
#' @export
release_path_sim <- function(x0 = 0, T = 1, alpha = 4, sigma = 1, theta = 0,
                        pts = 100, reps = 40, release_frac=.5){
  deltaT <- T/pts
  shiftT <- ceiling(release_frac * pts)
  ou_snap <- function(){
    x <- numeric(pts)
    x[1] <- x0
    for(i in 1:shiftT){
      mu <- x[i] * exp(-alpha * deltaT) + theta * (1 - exp(-alpha * deltaT))
      var <- sigma * sigma * (1 - exp(-2 * alpha * deltaT)) / (2 * alpha)
      x[i + 1] <- rnorm(1, mu, sqrt(var))
    }
    for(i in shiftT:pts){
      x[i + 1] <- rnorm(1, x[i], sigma*sqrt(deltaT))
    }
    x
  }
  X <- replicate(reps, ou_snap())
  out <- melt(X)
  names(out) <- c("time", "rep", "value")
  out
}


#' Release of constraint path simulator
#'
#' @param x0 starting position
#' @param T total time
#' @param sigma BM rate at start
#' @param sigma2 BM rate after shift
#' @param pts number of points per rep
#' @param reps number of replicates
#' @param release_frac fraction of time after which shift occurs
#' @return a melted data object ready for ggploting
#' @export
brownie_path_sim <- function(x0 = 0, T = 1, sigma = 1/sqrt(8), sigma2 = 1, 
                        pts = 100, reps = 40, release_frac=.5){
  deltaT <- T/pts
  shiftT <- ceiling(release_frac * pts)
  snap <- function(){
    x <- numeric(pts)
    x[1] <- x0
    for(i in 1:shiftT){
      x[i + 1] <- rnorm(1, x[i], sigma*sqrt(deltaT))
    }
    for(i in shiftT:pts){
      x[i + 1] <- rnorm(1, x[i], sigma2*sqrt(deltaT))
    }
    x
  }
  X <- replicate(reps, snap())
  out <- melt(X)
  names(out) <- c("time", "rep", "value")
  out
}


