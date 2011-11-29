# unit_test2.R
require(wrightscape)
data(labrids)

# Test of single OU
test <- dat[["prot.y"]]
test[!is.na(test)] <- rnorm(sum(!is.na(test)), sd=10)
ou  <- list(alpha = "global", sigma = "global", theta = "global")
m1 <- multiTypeOU(data=test, tree=tree, regimes=intramandibular, model_spec=ou)
test[!is.na(test)] <- rnorm(sum(!is.na(test)), sd=.1)
m2 <- multiTypeOU(data=test, tree=tree, regimes=intramandibular, model_spec=ou)
m1$alpha # larger variance has smaller alpha
m2$alpha


# test of the multi-case
testcase <- dat[["prot.y"]]
other.prot.data <- testcase[intramandibular=="other" & !is.na(testcase)]
intra.prot.data <- testcase[intramandibular!="other" & !is.na(testcase)] 
testcase[intramandibular=="other" & !is.na(testcase)] <- rnorm(length(other.prot.data), sd=.1)
testcase[intramandibular!="other" & !is.na(testcase)] <- rnorm(length(intra.prot.data), sd=10)
# Load internal functions
require(devtools)
load_all("..")
data <- testcase
lca <- lca_calc(tree)

# expect higher alpha on "other" (alpha2) to be best, higher on alpha1 is worst
# BUT NO!
multiOU_lik_lca(testcase, tree, intramandibular, alpha=c(2,2), sigma=c(5,5), theta=c(0.5, 0.5), Xo=0.5, lca)
multiOU_lik_lca(testcase, tree, intramandibular, alpha=c(5,2), sigma=c(5,5), theta=c(0.5, 0.5), Xo=0.5, lca)
multiOU_lik_lca(testcase, tree, intramandibular, alpha=c(2,5), sigma=c(5,5), theta=c(0.5, 0.5), Xo=0.5, lca)

a1_spec  <- list(alpha = "indep", sigma = "global", theta = "global")
a1 <- multiTypeOU(data=testcase, tree=tree, regimes=intramandibular, model_spec=a1_spec)
a1$alpha


