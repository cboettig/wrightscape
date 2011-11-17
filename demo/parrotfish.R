# mcmc_demo.R
rm(list=ls())
require(wrightscape)
require(snowfall)

data(parrotfish)

spec = list(alpha = "indep", sigma = "global", theta = "indep")
#traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")
#trait <- "prot.y"
traits <- c("close", "open", "gape.y",  "prot.y")

sfInit(par=T, 4)    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()

fits <- sfLapply(traits, function(trait){
  modelfit <- multiTypeOU(data=dat[trait], tree=tree, 
  regimes=intramandibular, model_spec=spec, control=list(maxit=2000)) 
  
  boots <- replicate(40, bootstrap(modelfit))
  out <- summary(modelfit, boots)
  out
})


## Plot ## 
regime.names <- levels(intramandibular)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) != 
     length(lower) | length(lower) != length(upper))
  stop("vectors must be same length")
  arrows(x,y + upper, x, y - lower, 
         angle = 90, code = 3, length=length, ...)
}

alphas <- sapply(fits, function(x)  x$Param.est["alpha",])
colnames(alphas) <- traits
alphas.se <- sapply(fits, function(x)   x$Param.SE["alpha",])
png("alphas.png", width=600)
  bars <- barplot(alphas, beside=T, main="alphas", legend.text=regime.names,
  ylim=c(0, max(alphas+alphas.se, na.rm=T)))
  error.bar(bars, alphas, alphas.se)
dev.off()

require(socialR)
upload("alphas.png", script="parrotfish.R", tag="phylogenetics")

