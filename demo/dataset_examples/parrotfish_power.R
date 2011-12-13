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
traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")

regimes <- intramandibular
  # declare function for shorthand
#sfInit(par=T, cpu=9)    # for debugging locally
#sfLibrary(wrightscape)
#sfLibrary(socialR)
#sfExportAll()

fits <- lapply(traits, function(trait){
	multi <- function(modelspec){ 
	 multiTypeOU(data = dat[[trait]], tree = tree, regimes = regimes, 
			    model_spec = modelspec, control = list(maxit=8000))

	}
	bm <- multi(list(alpha = "fixed", sigma = "indep", theta = "global")) 
	a1  <- multi(list(alpha = "indep", sigma = "global", theta = "global")) 
#	a2  <- multi(list(alpha = "indep", sigma = "global", theta = "indep")) 
#	full  <- multi(list(alpha = "indep", sigma = "indep", theta = "indep")) 
  mc <- montecarlotest(bm,a1)
  png(paste(trait, "_mc_parrotfish_", id, ".png", sep=""))
    plot(mc,show_data=TRUE, main=trait)
  dev.off()
#  upload("mc.png", gitaddr=gitaddr, tag="phylogenetics", comment=trait)
  mc
})

save(list=ls(), file=paste("parrotfish_power_", id, ".Rdat", sep=""))
