# unit_test3.R
rm(list=ls()) 
require(devtools)
load_all("..")
require(ouch)
data(bimac)
tree <- with(bimac, ouchtree(node, ancestor, time/max(time), species))
data <- log(bimac[['size']])
regimes <- bimac[["OU.LP"]]
combine <- (regimes=="small" | regimes=="large")
diff <-  combine & !is.na(data)
data[diff] <- rnorm(sum(diff), sd=5)  # Regime 1 should have the large variance
data[!combine & !is.na(data)] <- rnorm(sum(diff), mean=0, sd=.05)  # Regime 1 should have the large variance


regimes <- as.character(regimes)
regimes[combine] <- 1
regimes[!combine] <- 2
regimes <- as.factor(regimes)

lca <- lca_calc(tree)

theta <- c(mean(data[combine], na.rm=T), mean(data[!combine], na.rm=T))
Xo <- theta[match(regimes[1], levels(regimes))]  
# Here's the rub.  
bad <- multiOU_lik_lca(data, tree, regimes, alpha=c(4, 1e-7), sigma=12, theta=theta, Xo=Xo, lca=lca)
equal <- multiOU_lik_lca(data, tree, regimes, alpha=c(2,2), sigma=12, theta=theta, Xo=Xo, lca=lca)
good <- multiOU_lik_lca(data, tree, regimes, alpha=c(1e-7,4), sigma=12, theta=theta, Xo=Xo, lca=lca)
good > bad


## Show sigma's work
equal <- multiOU_lik_lca(data, tree, regimes, alpha=2, sigma=c(5,5), theta=3, Xo=3, lca)
bad <- multiOU_lik_lca(data, tree, regimes, alpha=2, sigma=c(5,10), theta=3, Xo=3, lca)
good <- multiOU_lik_lca(data, tree, regimes, alpha=2, sigma=c(10,5), theta=3, Xo=3, lca)
good > bad






## show alpha has correct effect globally
data[!is.na(data)] <- rnorm(sum(!is.na(data)), mean=3, sd=10)
bad <- multiOU_lik_lca(data, tree, regimes, alpha=4, sigma=10, theta=3, Xo=3, lca)
good <- multiOU_lik_lca(data, tree, regimes, alpha=.1, sigma=10, theta=3, Xo=3, lca)
good > bad


#ws <- wrightscape(data, tree, regimes, alpha=2.6, theta=3., sigma=.2, Xo=3.0)
#mt <- multiTypeOU(data, tree, regimes, alpha=2.6, theta=3., sigma=.2, Xo=3.0)

