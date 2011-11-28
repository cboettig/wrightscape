# File: labrid.R
# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD 

rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

# store the unique id of this script version
require(socialR)
gitaddr <- gitcommit("labrid.R")
id <- gitlog()$shortID


data(labrids)
traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")
#traits <- c("close", "open", "gape.y",  "prot.y")

sfInit(par=T, length(traits))    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()
fits <- sfLapply(traits, function(trait){
  alphas <- multiTypeOU(data=dat[trait], tree=tree, regimes=two_shifts, 
    model_spec = list(alpha = "indep", sigma = "global", theta = "indep"),
    alpha = c(.01, 10), sigma = c(.1, .1), control=list(maxit=1e5, tmax=20)) 
  alphas_out <- replicate(10, bootstrap(alphas))

  sigmas <- multiTypeOU(data=dat[trait], tree=tree, regimes=two_shifts, 
    model_spec = list(alpha = "global", sigma = "indep", theta = "indep"),
    sigma = c(.1, .1),
    control==list(maxit=1e5, tmax=20)) 
  sigmas_out <- replicate(10, bootstrap(sigmas))

  full <- multiTypeOU(data=dat[trait], tree=tree, regimes=two_shifts, 
    model_spec = list(alpha = "global", sigma = "indep", theta = "indep"), 
    alpha = c(.1, 10), sigma = c(.1, .1),
    control=list(maxit=2e5, tmax=20))   
    
  full_out <- replicate(10, bootstrap(full))

  list(ouma=alphas_out, oumv=sigmas_out, oumva=full_out)
})

names(fits) <- traits
data <- melt(fits)
names(data) <- c("regimes", "param", "rep", "value", "model", "trait")

p1 <- ggplot(subset(data,  param=="loglik")) + geom_boxplot(aes(model, value)) + facet_wrap(~ trait, scales="free_y")
p2 <- ggplot(subset(data, value < 100 & param %in% c("sigma", "alpha", "theta"))) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_grid(param ~ model, scales = "free") 
=======
fits <- lapply(traits, function(trait){

  # declare function for shorthand
  multi <- function(modelspec, reps = 20){
    m <- multiTypeOU(data = dat[trait], tree = tree, regimes = intramandibular, 
  		     model_spec = modelspec, 
#		     control = list(temp = 20, tmax = 50), method = "SANN"
		     control = list(maxit=5000)
		    ) 
    replicate(reps, bootstrap(m))
  }

  bm <- multi(list(alpha = "fixed", sigma = "indep", theta = "global"))
  s1 <- multi(list(alpha = "global", sigma = "indep", theta = "global")) 
  a1  <- multi(list(alpha = "indep", sigma = "global", theta = "global")) 
  s2 <- multi(list(alpha = "global", sigma = "indep", theta = "indep")) 
  a2  <- multi(list(alpha = "indep", sigma = "global", theta = "indep")) 

  list(bm=bm, s1=s1, a1=a1, s2=s2, a2=a2)
})

# Reformat and label data for plotting
names(fits) <- traits  # each fit is a different trait (so use it for a label)
data <- melt(fits)
names(data) <- c("regimes", "param", "rep", "value", "model", "trait")

#model likelihood
p1 <- ggplot(subset(data,  param=="loglik")) + 
      geom_boxplot(aes(model, value)) +
      facet_wrap(~ trait, scales="free_y")

# paramater estimates
p2 <- ggplot(subset(data, param %in% c("sigma", "alpha"))) +
      geom_boxplot(aes(trait, value, fill=regimes)) + 
      facet_grid(param ~ model, scales = "free") + scale_y_log() 

p3 <- ggplot(subset(data, param %in% c("sigma", "alpha"))) +
      geom_boxplot(aes(model, value, fill=regimes)) + 
      facet_grid(trait ~ param, scales = "free") 

save(list=ls(), file=sprintf("%s.Rdat", id))
ggsave(sprintf("%s_lik.png", id), p1)
ggsave(sprintf("%s_params_p2.png", id),  p2)
ggsave(sprintf("%s_params_p3.png", id),  p3)

save(list=ls(), file="parrotfish.Rdat")
require(ggplot2)
ggsave("parrotfish_lik.png", p1)
ggsave("parrotfish_params.png", p2)
require(socialR)
upload("parrotfish_*.png", script="parrotfish.R", tag="phylogenetics")

## For uploading plots at end  
#require(socialR); require(ggplot2); require(wrightscape)
#load(file=".Rdat") # must look up manually 
#upload(sprintf("%s_*.png", id), gitaddr=gitaddr, tag="phylogenetics")

#p <- ggplot(subset(data, param=="alpha" & value < 100)) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_grid(. ~ model) 
#ggplot(data) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_wrap(param~.)



