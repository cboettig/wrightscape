# mcmc_demo.R
rm(list=ls())
require(wrightscape)
require(snowfall)

source("labrid_data.R")
spec = list(alpha = "indep", sigma = "global", theta = "indep")
traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")

sfInit(par=T, cpu=4)
sfLibrary(wrightscape)
sfExportAll()

fits <- sfLapply(traits, function(trait){


  modelfit <- multiTypeOU(data=labrid$data[trait], tree=labrid$tree, 
  regimes=two_shifts, model_spec=spec) #,
# method ="SANN", control=list(maxit=100000,temp=50,tmax=20))


  reps <- sapply(1:100, function(i) simulate(modelfit))
  


  bootstrap <- sapply(1:40, 
    function(i){
      dat <- simulate(modelfit) 
      out <- update(modelfit, dat)
      names(out$alpha) <- paste("alpha", levels(modelfit$regimes), sep=".")
      names(out$sigma) <- paste("sigma", levels(modelfit$regimes), sep=".")
      names(out$theta) <- paste("theta", levels(modelfit$regimes), sep=".")
      pars <- rep(NA, 3*length(levels(modelfit$regimes)))
      if(out$convergence == 0) # only return values if successful
        pars <- c(out$alpha, out$sigma, out$theta)
      pars
    })

  est <- rbind(alpha = modelfit$alpha, sigma = modelfit$sigma,
               theta = modelfit$theta)
  SE <- sapply(1:dim(bootstrap)[1], function(i) sd(bootstrap[i,], na.rm=T) )
  SE <- t(matrix(SE, nrow = length(levels(modelfit$regimes))))
  rownames(SE) = c("alpha", "sigma", "theta") 
  colnames(SE) = levels(modelfit$regimes)
  colnames(est) = levels(modelfit$regimes)
  list(Param.est = est, Param.SE = SE)
})


regime.names <- levels(two_shifts)
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
#### Plot alphas ###
png("labrid_alphas.png", width=600)
  bars <- barplot(alphas, beside=T, main="alphas", legend.text=regime.names,
  ylim=c(0, max(alphas+alphas.se, na.rm=T)))
  error.bar(bars, alphas, alphas.se)
dev.off()

require(socialR)
upload("labrid_alphas.png", script="labrid_ou.R", tag="phylogenetics")

