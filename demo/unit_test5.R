#unit_test5.R

# Try on a star tree where we can use analytic results...
rm(list=ls())
require(wrightscape)
require(geiger)
data(labrids)
tree <- convert(lambdaTree(convert(tree), .01))
test <- dat[["prot.y"]]
regimes <- pharyngeal
levels(regimes) <- 1:length(levels(regimes))
wrasse <- !is.na(test) & (regimes == 1)
parrotfish <- !is.na(test) & (regimes == 2)
test[wrasse] <- rnorm(sum(wrasse), sd=1)
test[parrotfish] <- rnorm(sum(parrotfish), sd=10)
theta <- c(mean(test[regimes==1], na.rm=T), mean(test[regimes==2], na.rm=T))
Xo <- theta[match(regimes[1], levels(regimes))]  
# Should be approximately the optimum likelihood (along the alpha-sigma ridge)
sigma <- c(15,15)
alpha <- c(sigma[1]^2/(2*1), sigma[2]^2/(2*10))

require(devtools)
load_all("..")
multiOU_lik_lca(test, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo, lca=lca_calc(tree))



ou  <- list(alpha = "global", sigma = "global", theta = "global")
m1 <- multiTypeOU(data=test, tree=tree, regimes=regimes)



