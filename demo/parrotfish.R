# mcmc_demo.R
rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

data(parrotfish)

#traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")
traits <- c("close", "open")

sfInit(par=T, 2)    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()
fits <- sfLapply(traits, function(trait){
  alphas <- multiTypeOU(data=dat[trait], tree=tree, regimes=intramandibular, 
    model_spec = list(alpha = "indep", sigma = "global", theta = "indep"),
    alpha = c(.01, 10), sigma = c(.01, .01), control=list(maxit=5000)) 
  alphas_out <- replicate(20, bootstrap(alphas))

  sigmas <- multiTypeOU(data=dat[trait], tree=tree, regimes=intramandibular, 
    model_spec = list(alpha = "global", sigma = "indep", theta = "indep"), 
    control=list(maxit=5000))  
  sigmas_out <- replicate(20, bootstrap(sigmas))

  list(ouma=alphas_out, oumv=sigmas_out)
})
names(fits) <- traits
data <- melt(fits)
names(data) <- c("regimes", "param", "rep", "value", "model", "trait")


p <- ggplot(subset(data, value < 20 & param %in% c("sigma", "alpha", "theta"))) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_grid(param ~ model) 
print(p)

#p <- ggplot(subset(data, param=="alpha" & value < 100)) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_grid(. ~ model) 
#ggplot(data) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_wrap(param~.)




