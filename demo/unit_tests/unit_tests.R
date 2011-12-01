# Test case demonstrates sanity for wrightscape approach
rm(list=ls())
require(wrightscape)
data(labrids)
spec = list(alpha="indep", sigma="global", theta="global")
#spec = list(alpha="global", sigma="indep", theta="global")

testcase <- dat[["prot.y"]]
other.prot.data <- testcase[intramandibular=="other" & !is.na(testcase)]
intra.prot.data <- testcase[intramandibular!="other" & !is.na(testcase)] 

## Set the distribution of "other"
testcase[intramandibular=="other" & !is.na(testcase)] <- rnorm(length(other.prot.data), sd=.1)

#Set the distribution of the focal trait
testcase[intramandibular!="other" & !is.na(testcase)] <- rnorm(length(intra.prot.data), sd=10)


m <- multiTypeOU(data=testcase, tree=tree, regimes=intramandibular, model_spec=spec)
names(m$sigma) <- levels(intramandibular)
names(m$alpha) <- levels(intramandibular)

print("alphas")
print(m$alpha)

print("sigmas")
print(m$sigma)

#print("expected sd")
#print(sqrt(m$sigma^2 / 2*m$alpha))




rm(list=ls())
require(wrightscape)
data(labrids)
spec = list(alpha="indep", sigma="global", theta="global")
#spec = list(alpha="global", sigma="indep", theta="global")

testcase <- dat[["prot.y"]]
testcase[pharyngeal=="other" & !is.na(testcase)] -> other.prot.data
testcase[pharyngeal!="other" & !is.na(testcase)] -> intra.prot.data

## Set the distribution of "other"
testcase[pharyngeal=="other" & !is.na(testcase)] <- rnorm(length(other.prot.data), sd=.1)
#Set the distribution of the focal trait
testcase[pharyngeal!="other" & !is.na(testcase)] <- rnorm(length(intra.prot.data), sd=10)

m <- multiTypeOU(data=testcase, tree=tree, regimes=pharyngeal, model_spec=spec)
names(m$sigma) <- levels(pharyngeal)
names(m$alpha) <- levels(pharyngeal)

print("alphas")
print(m$alpha)
print("sigmas")
print(m$sigma)



