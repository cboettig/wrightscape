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
sfInit(par=T, cpu=4)    # for debugging locally
sfLibrary(wrightscape)
sfLibrary(socialR)
sfExportAll()

fits <- sfLapply(traits, function(trait){
	multi <- function(modelspec){ 
	 multiTypeOU(data = dat[["open"]], tree = tree, regimes = regimes, 
			    model_spec = modelspec, control = list(maxit=8000))

	}
#	bm <- multi(list(alpha = "fixed", sigma = "indep", theta = "global")) #uncensored
	bm2 <- multi(list(alpha = "fixed", sigma = "indep", theta = "indep")) #censored 
	a2  <- multi(list(alpha = "indep", sigma = "global", theta = "indep")) 
#	full  <- multi(list(alpha = "indep", sigma = "indep", theta = "indep")) 


  mc <- montecarlotest(bm2,a2)
  png("mc.png")
    plot(mc,show_data=TRUE)
  dev.off()

  upload("mc.png", gitaddr=gitaddr, tag="phylogenetics")
  mc
})

save(list=ls(), file="parrotfish_power.Rdat")
