
## CHECK branchlength times alpha
## CHECK lca_calc!!!


rm(list=ls())
require(wrightscape)
load("../data/labrids.rda")
regimes <- intramandibular


# Test of single OU
test <- dat[["prot.y"]]
test[!is.na(test)] <- rnorm(sum(!is.na(test)), sd=10)


require(devtools)
load_all("..")
lca <- lca_calc(tree)
sigma <- 10
alpha <- sigma^2/(2*10)
theta <- mean(test, na.rm=T)
Xo <- theta




ou  <- list(alpha = "global", sigma = "global", theta = "global")

m1 <- multiTypeOU(data=test, tree=tree, regimes=regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo, model_spec=ou)
tree <- convert(lambdaTree(convert(tree), .001))
m2 <- multiTypeOU(data=test, tree=tree, regimes=regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo, model_spec=ou)
print(m1$alpha) # larger variance has smaller alpha
print(m2$alpha)


