# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD 

rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

# store the unique id of this script version
require(socialR)
gitaddr <- gitcommit("parrotfish_power.R")
id <- gitlog()$shortID
print(id)

data(parrotfish)

regimes <- intramandibular
  # declare function for shorthand
sfInit(par=T, cpu=10)    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()

	multi <- function(modelspec){ 
	 multiTypeOU(data = dat[["open"]], tree = tree, regimes = regimes, 
			    model_spec = modelspec, control = list(maxit=8000))

	}
	s1 <- multi(list(alpha = "global", sigma = "indep", theta = "global")) 
	a1  <- multi(list(alpha = "indep", sigma = "global", theta = "global")) 
	s2 <- multi(list(alpha = "global", sigma = "indep", theta = "indep")) 
	a2  <- multi(list(alpha = "indep", sigma = "global", theta = "indep")) 


sfExportAll()
mc <- montecarlotest(s1,a1)
png("mc.png")
  plot(mc,show_data=TRUE)
dev.off()

upload("mc.png", gitaddr=gitaddr, tag="phylogenetics")
