#unit_test5.R
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

load_all("..")

## Equilibrium variance on non-star tree shouldn't be larger than on a star tree!!
test <- simulate_wrightscape(tree, regimes, alpha = 5, sigma=sqrt(10)*10, theta=0, Xo=0)$rep.1[[1]]
test[test == 0] <- NA
var(test, na.rm=T)
test <- simulate_wrightscape(star, regimes, alpha = 5, sigma=sqrt(10)*10, theta=0, Xo=0)$rep.1[[1]]
test[test == 0] <- NA
var(test, na.rm=T)





test <- dat[["prot.y"]]
wrasse <- !is.na(test) & (regimes == 1)
parrotfish <- !is.na(test) & (regimes == 2)
test[wrasse] <- rnorm(sum(wrasse), sd=1, mean=0)
test[parrotfish] <- rnorm(sum(parrotfish), sd=10, mean=0)


## IF YOU WANT A NON-STAR TREE, SHOULD SIMULATE TRAITS ON TREE 
test <- simulate_wrightscape(tree, regimes, alpha = c(5, .5), sigma=c(10,10), theta=c(0,0), Xo=0)$rep.1[[1]]
test[test == 0] <- NA
var1 <- var(test[wrasse], na.rm=T)
var2 <- var(test[parrotfish], na.rm=T)
print(c(var1, var2))


theta <- c(mean(test[regimes==1], na.rm=T), mean(test[regimes==2], na.rm=T))
Xo <- theta[match(regimes[1], levels(regimes))]  
# Should be approximately the optimum likelihood (along the alpha-sigma ridge)
sigma <- c(20,20)
alpha <- c(sigma[1]^2/(2*var1), sigma[2]^2/(2*var2))

require(devtools)
load_all("..")
good <- multiOU_lik_lca(test, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo, lca=lca_calc(tree), scale=max(tree@times)*0.01)
bad <- multiOU_lik_lca(test, tree, regimes, alpha=rev(alpha), sigma=sigma, theta=theta, Xo=Xo, lca=lca_calc(tree), scale=max(tree@times)*0.01)
good
bad
print(good > bad)


a1 <- list(alphas="indep", sigmas="global", thetas="global")
a2 <- list(alphas="indep", sigmas="global", thetas="indep")
m1 <- multiTypeOU(test, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo, model_spec=a2)
print(m1$alpha) # YAY THIS IS CORRECT!!

m2 <- multiTypeOU(test, tree, regimes, model_spec=a2)
print(m2$alpha)


m3 <-  multiTypeOU(test, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo)
print(m3$alpha)


