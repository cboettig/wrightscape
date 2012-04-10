


require(wrightscape)
require(ggplot2)

plot.path_sim <- function(x, ...)
  ggplot(x) + geom_line(aes(time, value, group = rep), alpha = .05)


  reps <- 500
  X <- bm_path_sim(reps = reps)
  Y <- ou_path_sim(reps = reps, alpha = 6)
  Z <- release_path_sim(reps = reps, alpha = 6, release_frac=.7)
  W <- brownie_path_sim(reps = reps, sigma = sqrt(1/.7)/sqrt(2*6), sigma2=1, release_frac=.7 ) # has 1/release_frac time = 2

  px <- plot.path_sim(X) + opts(title = "Brownian Motion") + coord_cartesian(ylim = c(-2,2))
  py <- plot.path_sim(Y) + opts(title = "Ornstein-Uhlenbeck") + coord_cartesian(ylim = c(-2,2))
  pz <- plot.path_sim(Z) + opts(title = "Release of Constraint") + coord_cartesian(ylim = c(-2,2))
  pw <- plot.path_sim(W) + opts(title = "Accelerated Evolution") + coord_cartesian(ylim = c(-2,2))

#  grid.newpage()
  cairo_pdf("release.pdf")
  pushViewport(viewport(layout = grid.layout(2,2)))
  vplayout <- function(x, y) viewport(layout.pos.row = x,
  layout.pos.col = y)
  print(px, vp = vplayout(1, 1))
  print(py, vp = vplayout(1, 2))
  print(pw, vp = vplayout(2, 1))
  print(pz, vp = vplayout(2, 2))
  dev.off()




