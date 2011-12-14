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

#traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")
traits <- c("close", "open", "kt", "gape.y", "AM.y")
regimes <- two_shifts 

  # declare function for shorthand
sfInit(par=T, cpu=5)    # for debugging locally
fits <- lapply(traits, function(trait){
	multi <- function(modelspec){ 
	 multiTypeOU(data = dat[[trait]], tree = tree, regimes = regimes, 
			    model_spec = modelspec, control = list(maxit=8000))

	}
	bm2 <- multi(list(alpha = "fixed", sigma = "indep", theta = "global")) 
	a2  <- multi(list(alpha = "indep", sigma = "global", theta = "global")) 
	a2t2  <- multi(list(alpha = "indep", sigma = "global", theta = "indep")) 
	t2  <- multi(list(alpha = "global", sigma = "global", theta = "indep")) 

  mc <- montecarlotest(bm2,a2)
  bm_a2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(bm,a2t2)
  bm_a2t2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(a2,a2t2)
  a2_a2t2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(t2,a2)
  t2_a2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(t2,a2t2)
  t2_a2t2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))

  list(bm_a2=bm_a2, bm_a2t2=bm_a2t2, a2_a2t2=a2_a2t2, t2_a2=t2_a2, t2_a2t2=t2_a2t2)
})


names(fits) <- traits
dat <- melt(fits)

names(dat) <- c("value", "type", "comparison", "trait")

require(ggplot2)

p1 <- ggplot(subset(dat, abs(value) < 1e3)) + 
      geom_boxplot(aes(type, value)) +
      facet_grid(comparison ~ trait, scales="free")
ggsave("modelchoice.png", p1)


save(list=ls(), file=paste("labrid_power_", id, ".Rdat", sep=""))
