#unit_test5.R

# Try on a star tree where we can use analytic results...
require(wrightscape)
require(geiger)
#data(labrids)
load("../data/labrids.rda")
tree <- convert(lambdaTree(convert(tree), .1))
test <- dat[["prot.y"]]
regimes <- pharyngeal
levels(regimes) <- 1:length(levels(regimes))
wrasse <- !is.na(test) & (regimes == 1)
parrotfish <- !is.na(test) & (regimes == 2)
test[wrasse] <- rnorm(sum(wrasse), sd=1, mean=0)
test[parrotfish] <- rnorm(sum(parrotfish), sd=10, mean=0)
theta <- c(mean(test[regimes==1], na.rm=T), mean(test[regimes==2], na.rm=T))
Xo <- theta[match(regimes[1], levels(regimes))]  
# Should be approximately the optimum likelihood (along the alpha-sigma ridge)
sigma <- c(10,10)
alpha <- c(sigma[1]^2/(2*1), sigma[2]^2/(2*10))

require(devtools)
load_all("..")
good <- multiOU_lik_lca(test, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo, lca=lca_calc(tree))
bad <- multiOU_lik_lca(test, tree, regimes, alpha=rev(alpha), sigma=sigma, theta=theta, Xo=Xo, lca=lca_calc(tree))
print(good > bad)


a1 <- list(alphas="indep", sigmas="global", thetas="global")
a2 <- list(alphas="indep", sigmas="global", thetas="indep")
m1 <- multiTypeOU(test, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo, model_spec=a2)
print(m1$alpha) # YAY THIS IS CORRECT!!

m2 <- multiTypeOU(test, tree, regimes, model_spec=a2)
print(m2$alpha)


m3 <-  multiTypeOU(test, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo)
print(m3$alpha)


