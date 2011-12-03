# File: labrid.R
# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD 

rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

# store the unique id of this script version
require(socialR)
gitaddr <- gitcommit("labrids.R")
id <- gitlog()$shortID

print(id)


data(labrids)


  # declare function for shorthand
multi <- function(modelspec) 
	multiTypeOU(data = dat[["open"]], tree = tree, regimes = intramandibular, 
		    model_spec = modelspec, control = list(maxit=3000)) 

s2 <- multi(list(alpha = "global", sigma = "indep", theta = "indep")) 
a2  <- multi(list(alpha = "indep", sigma = "global", theta = "indep")) 


sfInit(par=T, 10)    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()
out <- montecarlotest(s2,a2,nboot=100)
