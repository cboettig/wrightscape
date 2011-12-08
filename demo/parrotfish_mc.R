# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD 

rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

# store the unique id of this script version
require(socialR)
gitaddr <- gitcommit("parrotfish_mc.R")
id <- gitlog()$shortID

print(id)


data(parrotfish)
traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")

regimes <- intramandibular
  # declare function for shorthand
sfInit(par=T, 4)    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()

fits <- sfLapply(traits, function(trait){
	multi <- function(modelspec){ 
		out <- multiTypeOU(data = dat[[trait]], tree = tree, regimes = regimes, 
			    model_spec = modelspec, control = list(maxit=8000))

	      ## plots require this formatting.  should make as a seperate function
	      ## so can return full model fits.   
	      n <- length(levels(out$regimes))
	      Xo <- rep(out$Xo,n) 
	      loglik <- rep(out$loglik, n)
	      pars <- cbind(out$alpha, out$sigma, out$theta, Xo, loglik)
	      rownames(pars) <- levels(out$regimes)
	      colnames(pars) <- c("alpha", "sigma", "theta", "Xo", "loglik")
	      if(out$convergence != 0) # only return values if successful
		pars[,] <- NA
	      pars
	}

  # no point in bm with theta indep, since the only theta that matters is root value
	bm <- multi(list(alpha = "fixed", sigma = "indep", theta = "global")) 
	a1  <- multi(list(alpha = "indep", sigma = "global", theta = "indep")) 
	a2  <- multi(list(alpha = "indep", sigma = "global", theta = "indep")) 
	full  <- multi(list(alpha = "indep", sigma = "indep", theta = "indep")) 
	list(bm=bm,a1=a1,a2=a2,full=full)
})

# Reformat and label data for plotting
names(fits) <- traits  # each fit is a different trait (so use it for a label)
data <- melt(fits)
names(data) <- c("regimes", "param", "value", "model", "trait")


#model likelihood
p1 <- ggplot(subset(data,  param=="loglik")) + geom_boxplot(aes(model, value)) +
      facet_wrap(~ trait, scales="free_y")
p2 <-  ggplot(subset(data, param %in% c("alpha") & model %in% c("a1", "a2", "full") ),
              aes(model, value, fill=regimes)) + geom_bar(position="dodge") +  
       facet_wrap(~trait, scales="free_y")
ggsave(sprintf("%s_p1.png", id), p1)
ggsave(sprintf("%s_p2.png", id),  p2)
require(socialR)
upload(sprintf("%s_p*.png", id), gitaddr=gitaddr, tag="phylogenetics")

save(list=ls(), file=sprintf("%s.Rdat", id))
print(id)

