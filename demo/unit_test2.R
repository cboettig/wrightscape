# unit_test2.R
rm(list=ls())
require(wrightscape)
data(labrids)
regimes <- intramandibular

# test of the multi-case
test <- dat[["prot.y"]]
group_1 <- (regimes==levels(regimes)[1]) & !is.na(test)
group_2 <- (regimes==levels(regimes)[2]) & !is.na(test)
test[group_1] <- rnorm(sum(group_1), sd=10)
test[group_2] <- rnorm(sum(group_2), sd=1)

# Load internal functions
require(devtools)
load_all("..")
lca <- lca_calc(tree)

sigma <- c(10,10)
alpha <- c(sigma[1]^2/(2*1), sigma[2]^2/(2*10))

theta <- c(mean(test[regimes==levels(regimes)[1]], na.rm=T), 
	   mean(test[regimes==levels(regimes)[2]], na.rm=T))
Xo <- theta[match(regimes[1], levels(regimes))]  

# expect higher alpha on "other" (alpha2) to be best, higher on alpha1 is worst
even <- multiOU_lik_lca(test, tree, regimes,  alpha=mean(alpha), sigma=sigma, theta=theta, Xo=Xo, lca=lca_calc(tree))
bad <- multiOU_lik_lca(test, tree, regimes, alpha=c(alpha[2], alpha[1]), sigma=sigma, theta=theta, Xo=Xo, lca=lca_calc(tree))
good <- multiOU_lik_lca(test, tree, regimes, alpha=alpha, sigma=sigma, theta=theta, Xo=Xo, lca=lca_calc(tree) )

good > bad
# BUT NO!  FIXME -- check out means of the two groups!!


multiOU_lik_lca(testcase, tree, intramandibular, alpha=c(2,2), sigma=c(5,5), theta=c(0.5, 0.5), Xo=0.5, lca)
multiOU_lik_lca(testcase, tree, intramandibular, alpha=c(5,2), sigma=c(5,5), theta=c(0.5, 0.5), Xo=0.5, lca)
multiOU_lik_lca(testcase, tree, intramandibular, alpha=c(2,5), sigma=c(5,5), theta=c(0.5, 0.5), Xo=0.5, lca)

## MODEL FIT
a1_spec  <- list(alpha = "indep", sigma = "global", theta = "global")
a1 <- multiTypeOU(data=testcase, tree=tree, regimes=intramandibular, model_spec=a1_spec)
a1$alpha





