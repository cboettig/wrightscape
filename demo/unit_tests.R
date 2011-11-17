# Test case demonstrates sanity for wrightscape approach
rm(list=ls())
require(wrightscape)
source("parrotfish_data.R")
spec = list(alpha="fixed", sigma="indep", theta="global")

testcase <- labrid$data[["prot.y"]]
testcase[intramandibular=="other" & !is.na(testcase)] -> other.prot.data
testcase[intramandibular!="other" & !is.na(testcase)] -> intra.prot.data
testcase[intramandibular=="other" & !is.na(testcase)] <- rnorm(length(other.prot.data), sd=1)
testcase[intramandibular!="other" & !is.na(testcase)] <- rnorm(length(intra.prot.data), sd=.1)
modelfit <- multiTypeOU(data=testcase, tree=labrid$tree, regimes=intramandibular, model_spec=spec)

names(modelfit$sigma) <- levels(intramandibular)
modelfit$sigma
