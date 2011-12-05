# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD 

rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

# store the unique id of this script version
require(socialR)
gitaddr <- gitcommit("labrid_power.R")
id <- gitlog()$shortID
print(id)

data(labrids)

regimes <- pharyngeal
  # declare function for shorthand
sfInit(par=T, 10)    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()

	multi <- function(modelspec){ 
	 multiTypeOU(data = dat[["close"]], tree = tree, regimes = regimes, 
			    model_spec = modelspec, control = list(maxit=8000))

	}
	s1 <- multi(list(alpha = "global", sigma = "indep", theta = "global")) 
	a1  <- multi(list(alpha = "indep", sigma = "global", theta = "global")) 

sfExportAll()
mc <- montecarlotest(s1,a1)


