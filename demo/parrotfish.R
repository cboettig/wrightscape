# mcmc_demo.R
rm(list=ls())
require(wrightscape)
require(snowfall)

data(parrotfish)

spec = list(alpha = "global", sigma = "indep", theta = "indep")
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



# value   variable   regime   model
# 0.011   alpha.e    intra    oumv
# 0.004   alpha.se   intra    oumv
# 1.124   sigma      intra    oumv 

# rep   value   variable   regime   model
# 1     0.011   alpha      intra    oumv
# 2     0.014   alpha      intra    oumv
# 1     1.124   sigma      intra    oumv 





param <- "sigma"
regime.names <- levels(intramandibular)


# a handy function for errorbars on bar plots
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) != 
     length(lower) | length(lower) != length(upper))
  stop("vectors must be same length")
  arrows(x,y + upper, x, y - lower, 
         angle = 90, code = 3, length=length, ...)
}

# collate the data 
estimate <- sapply(fits, function(x)  x$Param.est[param,])
estimate.se <- sapply(fits, function(x)   x$Param.SE[param,])
colnames(estimate) <- traits
colnames(estimate.se) <- traits

print(estimate)
print(estimate.se)

## Plot ## 
png("param.png", width=600)
  bars <- barplot(estimate, beside=T, main=param, legend.text=regime.names,
  ylim=c(0, max(estimate+estimate.se, na.rm=T)))
  error.bar(bars, estimate, estimate.se)
dev.off()

## share 
require(socialR)
upload("param.png", script="parrotfish.R", tag="phylogenetics")

