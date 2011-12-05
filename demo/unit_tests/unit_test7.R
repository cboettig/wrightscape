#unit_test7.R
rm(list=ls())
# Try on a star tree where we can use analytic results...
require(wrightscape)
require(geiger)
#data(labrids)
load("../data/labrids.rda")
star <- convert(lambdaTree(convert(tree), .001))
#tree <- tr
regimes <- pharyngeal
levels(regimes) <- 1:length(levels(regimes))

require(devtools)
load_all("..")

## Equilibrium variance on non-star tree shouldn't be larger than on a star tree!!
test <- simulate_wrightscape(tree, regimes, alpha = 5, sigma=sqrt(10)*10, theta=0, Xo=0)$rep.1[[1]]
test[test == 0] <- NA
print( var(test, na.rm=T) )
test <- simulate_wrightscape(star, regimes, alpha = 5, sigma=sqrt(10)*10, theta=0, Xo=0)$rep.1[[1]]
test[test == 0] <- NA
print( var(test, na.rm=T) )


