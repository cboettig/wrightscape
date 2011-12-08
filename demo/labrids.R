# File: labrid.R
# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD 

rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

# store the unique id of this script version
require(socialR)
gitaddr <- gitcommit("labrids.R")
id <- gitlog()$shortID

print(id)


data(labrids)
#traits <- c("bodymass", "close", "open", "kt", "gape.x",  "prot.x", "AM.x", "SH.x", "LP.x")
traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")

regimes <- intramandibular

sfInit(par=T, 9)    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()
fits <- sfLapply(traits, function(trait){

  # declare function for shorthand
  multi <- function(modelspec, reps = 20){
    m <- multiTypeOU(data = dat[trait], tree = tree, regimes = regimes, 
  		     model_spec = modelspec, 
		     control = list(maxit=3000)
		    ) 
    replicate(reps, bootstrap(m))
  }

  bm <- multi(list(alpha = "fixed", sigma = "indep", theta = "global"))
  a1  <- multi(list(alpha = "indep", sigma = "global", theta = "global")) 
  a2  <- multi(list(alpha = "indep", sigma = "global", theta = "indep")) 

  list(bm=bm, a1=a1, a2=a2)
})

# Reformat and label data for plotting
names(fits) <- traits  # each fit is a different trait (so use it for a label)
data <- melt(fits)
names(data) <- c("regimes", "param", "rep", "value", "model", "trait")

#model likelihood
p1 <- ggplot(subset(data,  param=="loglik")) + 
      geom_boxplot(aes(model, value)) +
      facet_wrap(~ trait, scales="free_y")

p2 <-  ggplot(subset(data, param %in% c("sigma", "alpha") & model != "bm"),
              aes(model, value, fill=regimes)) + 
#       stat_summary(fun.y=mean, geom="bar", position="dodge") + # add bars for some extra ink...
       stat_summary(fun.data=median_hilow, geom="pointrange", aes(color=regimes), 
                    position = position_dodge(width=0.90), conf.int=.5) +
       scale_y_log() + 
       facet_grid(param ~ trait, scales = "free_y")  
#       coord_cartesian(ylim=c(0,10), wise=TRUE)


## Just plot the parameters seperately
p3 <-  ggplot(subset(data, param %in% c("sigma") ), 
              aes(model, value, fill=regimes)) + 
       stat_summary(fun.y=mean, geom="bar", position="dodge", alpha=.5) + # add bars for some extra ink...
       stat_summary(fun.data=median_hilow, geom="pointrange", aes(color=regimes), 
                    position = position_dodge(width=0.90), conf.int=.5) +
#       scale_y_log() + 
       facet_grid(param ~ trait, scales = "free_y") + 
#       coord_cartesian(ylim=c(0,4), wise=TRUE) +  # easiest to just adjust the zoom limits manually still...
       opts(title="sigma")

p4 <-  ggplot(subset(data, param %in% c("alpha") ), 
              aes(model, value, fill=regimes)) + 
       stat_summary(fun.y=mean, geom="bar", position="dodge", alpha=.5) + # add bars for some extra ink...
       stat_summary(fun.data=median_hilow, geom="pointrange", aes(color=regimes), 
                    position = position_dodge(width=0.90), conf.int=.5) +
#       scale_y_log() + 
       facet_grid(param ~ trait, scales = "free_y") + 
#       coord_cartesian(ylim=c(0,10), wise=TRUE) +
        opts(title = "alpha")



save(list=ls(), file=sprintf("%s.Rdat", id))
ggsave(sprintf("%s_lik.png", id), p1)
ggsave(sprintf("%s_params_p2.png", id),  p2)
ggsave(sprintf("%s_params_p3.png", id),  p3)
ggsave(sprintf("%s_params_p4.png", id),  p4)


print(id)

## For uploading plots at end  
#require(socialR); require(ggplot2); require(wrightscape)
#load(file=".Rdat") # must look up manually 
#upload(sprintf("%s_*.png", id), gitaddr=gitaddr, tag="phylogenetics")

