
An example using ensemble simulations to visualize the difference between the models.  

``` {r }
require(wrightscape)
require(ggplot2)
require(reshape2)
require(grid)
````

Plotting function
``` {r }
plot.path_sim <- function(x, ...)
  ggplot(x) + geom_line(aes(time, value, group = rep), alpha = .05)
````

Simulations
``` {r }
reps <- 500
X <- bm_path_sim(reps = reps)
Y <- ou_path_sim(reps = reps, alpha = 6)
Z <- release_path_sim(reps = reps, alpha = 6, release_frac=.7)
W <- brownie_path_sim(reps = reps, sigma = sqrt(1/.7)/sqrt(2*6), sigma2=1, release_frac=.7 ) # has 1/release_frac time = 2
````

Create the plots 

``` {r fig.keep="last", fig.show="hold"}
px <- plot.path_sim(X) + opts(title = "Brownian Motion") + coord_cartesian(ylim = c(-2,2))
py <- plot.path_sim(Y) + opts(title = "Ornstein-Uhlenbeck") + coord_cartesian(ylim = c(-2,2))
pz <- plot.path_sim(Z) + opts(title = "Release of Constraint") + coord_cartesian(ylim = c(-2,2))
pw <- plot.path_sim(W) + opts(title = "Accelerated Evolution") + coord_cartesian(ylim = c(-2,2))
pushViewport(viewport(layout = grid.layout(2,2)))
vplayout <- function(x, y) viewport(layout.pos.row = x,
layout.pos.col = y)
print(px, vp = vplayout(1, 1))
print(py, vp = vplayout(1, 2))
print(pw, vp = vplayout(2, 1))
print(pz, vp = vplayout(2, 2))
````



