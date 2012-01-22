#' brownian motion path simulator
#'
#'
bm_path_sim <- function(x0 = 0, T = 1, sigma = 1, pts = 100, reps = 40){
  deltaT <- T/pts
  X <- replicate(reps,  x0 + cumsum(rnorm(pts, 0, sigma*sqrt(deltaT))))
  out <- melt(X)
  names(out) <- c("time", "rep", "value")
  out
}


#' OU path simulator
ou_path_sim <- function(x0 = 0, T = 1, alpha = 2, sigma = 1, theta = 0,
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
release_path_sim <- function(x0 = 0, T = 1, alpha = 2, sigma = 1, theta = 0,
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



plot.path_sim <- function(x, ...)
  ggplot(x) + geom_line(aes(time, value, group = rep), alpha = .3)


example_path <- function(){
  require(ggplot2)
  reps <- 100
  X <- bm_path_sim(reps = reps)
  Y <- ou_path_sim(reps = reps, alpha = 4)
  Z <- release_path_sim(reps = reps, alpha = 4)

  px <- plot.path_sim(X) + opts(title = "Brownian Motion") + coord_cartesian(ylim = c(-2,2))
  py <- plot.path_sim(Y) + opts(title = "Ornstein-Uhlenbeck") + coord_cartesian(ylim = c(-2,2))
  pz <- plot.path_sim(Z) + opts(title = "Release of Constraint") + coord_cartesian(ylim = c(-2,2))


  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,3)))
  vplayout <- function(x, y) viewport(layout.pos.row = x,
  layout.pos.col = y)
  print(px, vp = vplayout(1, 1))
  print(pz, vp = vplayout(1, 2))
  print(py, vp = vplayout(1, 3))
}

