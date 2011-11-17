# mcmc_demo.R
rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

data(parrotfish)

#traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")
traits <- c("close", "open", "gape.y",  "prot.y")

sfInit(par=T, 4)    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()
fits <- sfLapply(traits, function(trait){
  alphas <- multiTypeOU(data=dat[trait], tree=tree, regimes=intramandibular, 
    model_spec = list(alpha = "indep", sigma = "global", theta = "indep"),
    alpha = c(.01, 10), sigma = c(.01, .01), control=list(maxit=1e5, tmax=20), method="SANN") 
  alphas_out <- replicate(10, bootstrap(alphas))

  sigmas <- multiTypeOU(data=dat[trait], tree=tree, regimes=intramandibular, 
    model_spec = list(alpha = "global", sigma = "indep", theta = "indep"), 
    control==list(maxit=1e5, tmax=20), method="SANN") 
  sigmas_out <- replicate(10, bootstrap(sigmas))

  full <- multiTypeOU(data=dat[trait], tree=tree, regimes=intramandibular, 
    model_spec = list(alpha = "global", sigma = "indep", theta = "indep"), 
    control=list(maxit=2e5, tmax=20), method="SANN")   
  full_out <- replicate(10, bootstrap(full))

  list(ouma=alphas_out, oumv=sigmas_out, oumva=full_out)
})

names(fits) <- traits
data <- melt(fits)
names(data) <- c("regimes", "param", "rep", "value", "model", "trait")

p1 <- ggplot(subset(data,  param=="loglik")) + geom_boxplot(aes(model, value)) + facet_wrap(~ trait, scales="free_y")
p2 <- ggplot(subset(data, value < 100 & param %in% c("sigma", "alpha", "theta"))) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_grid(param ~ model, scales = "free") 

save(list=ls(), file="parrotfish.Rdat")
require(ggplot2)
ggsave("parrotfish_lik.png", p1)
ggsave("parrotfish_params.png", p2)
require(socialR)
upload("parrotfish_*.png", script="parrotfish.R", tag="phylogenetics")


#p <- ggplot(subset(data, param=="alpha" & value < 100)) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_grid(. ~ model) 
#ggplot(data) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_wrap(param~.)




