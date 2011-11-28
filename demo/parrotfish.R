# mcmc_demo.R
rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

# store the unique id of this script version
require(socialR)
gitaddr <- gitcommit("parrotfish.R")
id <- gitlog()$shortID



data(parrotfish)
#traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")
traits <- c("close", "open", "gape.y",  "prot.y")

sfInit(par=T, 4)    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()
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

p1 <- ggplot(subset(data,  param=="loglik")) + 
      geom_boxplot(aes(model, value)) +
      facet_wrap(~ trait, scales="free_y")

p2 <- ggplot(subset(data, value < 100 & param %in% c("sigma", "alpha", "theta"))) +
      geom_boxplot(aes(trait, value, fill=regimes)) + 
      facet_grid(param ~ model, scales = "free") 

save(list=ls(), file=sprintf("%s.Rdat", id))
ggsave(sprintf("%s_lik.png", id), p1)
ggsave(sprintf("%s_params.png", id),  p2)


require(socialR)
#load(file=sprintf("%s.Rdat", id))
#upload(sprintf("%s_*.png", id), gitaddr=gitaddr, tag="phylogenetics")


#p <- ggplot(subset(data, param=="alpha" & value < 100)) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_grid(. ~ model) 
#ggplot(data) + geom_boxplot(aes(trait, value, fill=regimes)) + facet_wrap(param~.)




