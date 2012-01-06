# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD 

rm(list=ls())
require(wrightscape)
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
fits <- lapply(traits, function(trait){
	multi <- function(modelspec){ 
	 multiTypeOU(data = dat[[trait]], tree = tree, regimes = regimes, 
			    model_spec = modelspec, control = list(maxit=8000))

	}
	bm <- multi(list(alpha = "fixed", sigma = "global", theta = "global")) 
	ou <- multi(list(alpha = "global", sigma = "global", theta = "global")) 
	bm2 <- multi(list(alpha = "fixed", sigma = "indep", theta = "global")) 
	a2  <- multi(list(alpha = "indep", sigma = "global", theta = "global")) 
	t2  <- multi(list(alpha = "global", sigma = "global", theta = "indep")) 

  mc <- montecarlotest(bm2,a2)
  bm2_a2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(bm,ou)
  bm_ou <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(bm,bm2)
  bm_bm2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(ou,bm2)
  ou_bm2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(t2,a2)
  t2_a2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))
  mc <- montecarlotest(bm2,t2)
  bm2_t2 <- list(null=mc$null_dist, test=mc$test_dist, lr=-2*(mc$null$loglik-mc$test$loglik))

  list(brownie_vs_alphas=bm2_a2, brownie_vs_thetas=bm2_t2, thetas_vs_alphas=t2_a2,
       bm_vs_brownie=bm_bm2,  bm_vs_ou=bm_ou, ou_vs_brownie=ou_bm2)
})

names(fits) <- traits
dat <- melt(fits)
names(dat) <- c("value", "type", "comparison", "trait")

save(list=ls(), file=paste("labrid_power_", id, ".Rdat", sep=""))


require(ggplot2)
r <- cast(dat, comparison ~ trait, function(x) quantile(x, c(.10,.90)))
subdat <- subset(dat, abs(value) < max(abs(as.matrix(r))))

p1 <- ggplot(subdat) + 
      geom_boxplot(aes(type, value)) +
      facet_grid(trait ~ comparison, scales="free_y") 
ggsave(paste(id, "_modelchoice.png", sep=""), p1)

## Tough to see everything on such a grid
for(tr in traits){
  p <- ggplot(subset(subdat, trait==tr)) +  geom_boxplot(aes(type, value)) +   facet_wrap(~ comparison, scales="free_y")
  ggsave(paste(id, "_", tr, ".png", sep=""), p)
}


