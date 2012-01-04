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

traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")
#traits <- c("close", "open", "kt", "gape.y", "AM.y")
regimes <- two_shifts 

  # declare function for shorthand
sfInit(par=F)    # for debugging locally
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
  mc <- montecarlotest(bm2,a2t2)
  bm_a2t2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(a2,a2t2)
  a2_a2t2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(t2,a2)
  t2_a2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(t2,a2t2)
  t2_a2t2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(bm2,t2)
  bm_t2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))

  list(sigmas_vs_alphas=bm_a2, sigmas_vs_thetas=bm_t2, thetas_vs_alphas=t2_a2,
       sigmas_vs_alphasthetas=bm_a2t2,  alphas_vs_alphasthetas=a2_a2t2,  
       thetas_vs_alphasthetas=t2_a2t2)
})

names(fits) <- traits
dat <- melt(fits)
names(dat) <- c("value", "type", "comparison", "trait")

save(list=ls(), file=paste("labrid_power_", id, ".Rdat", sep=""))


require(ggplot2)
#r <- cast(dat, comparison ~ trait, function(x) quantile(x, c(.05,.95)))
#subdat <- subset(dat, abs(value) < max(abs(as.matrix(r))))

p1 <- ggplot(subdat) + 
      geom_boxplot(aes(type, value)) +
      facet_grid(trait ~ comparison, scales="free_y") 
ggsave(paste(id, "_modelchoice.png", sep=""), p1)

## Tough to see everything on such a grid
for(tr in traits){
  p <- ggplot(subset(dat, trait==tr)) +  geom_boxplot(aes(type, value)) +   facet_wrap(~ comparison, scales="free_y")
  ggsave(paste(id, "_", tr, ".png", sep=""), p)
}


