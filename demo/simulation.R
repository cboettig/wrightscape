# simulation.R
rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

# store the unique id of this script version
require(socialR)
gitaddr <- gitcommit("simulation.R")
id <- gitlog()$shortID





data(labrids)
#
#  Create a dataset by simulation, where the parrotfish all have lower alpha
#
a1_spec  <- list(alpha = "indep", sigma = "global", theta = "global")
a1 <- multiTypeOU(data = dat["close"], tree = tree, regimes = pharyngeal, 
	     model_spec = a1_spec,  control = list(maxit=5000))
names(a1$alpha) <- levels(pharyngeal) 
a1$alpha["other"] <- 10
a1$alpha["pharyngeal"] <- 1e-10
#a1$sigma <- c(15, 15)  # We can keep those parameters estimated from data or update them
#a1$theta <- c(0,0)   
dat[["simulated_a1"]] <-simulate(a1)[[1]]

# Check out the variance in the relative groups -- it should be larger in parrotfish
# in an extreme example, but need not be... 
regime <- pharyngeal
testcase <- dat[["simulated_a1"]]
testcase[regime =="other" & !is.na(testcase) & testcase != 0] -> lowvar
testcase[regimei !="other" & !is.na(testcase) & testcase != 0] -> highvar
var(lowvar)
var(highvar)

# estimate a model from this simulated data
dummy <- update(a1, dat[["simulated_a1"]]) 
print(dummy$alpha) # let's see what kind of 




## We can repeat the whole thing with a model based on differnt sigmas, to make sure 

s1_spec  <- list(alpha = "global", sigma = "indep", theta = "global")
s1 <- multiTypeOU(data = dat["close"], tree = tree, regimes = pharyngeal, 
	     model_spec = s1_spec,  control = list(maxit=5000))
# Order of entries in sigma is the order regime names are given by levels (alphabetical)
names(s1$sigma) <- levels(pharyngeal)
s1$sigma["other"] <- .1
s1$sigma["pharyngeal"] <- 10
testcase <- simulate(s1)[[1]]
testcase[pharyngeal=="other" & !is.na(testcase) & testcase != 0 ] -> lowvar
testcase[pharyngeal!="other" & !is.na(testcase) & testcase != 0 ] -> highvar
var(lowvar)
var(highvar)
# and the recovered estimates are pretty close to the underlying simulation parameters
dummy <- update(s1, testcase)
dummy$sigma
dat[["simulated_s1"]] <-simulate(s1)[[1]]


## Now we have a trait where change in alpha is responsible, 
## and one in which sigma change is responsible. 
## Can we correctly identify each??

traits <- c("simulated_a1", "simulated_s1")
sfInit(par=T, 2)    # for debugging locally
sfLibrary(wrightscape)
sfExportAll()
fits <- lapply(traits, function(trait){
  # declare function for shorthand
  multi <- function(modelspec, reps = 20){
    m <- multiTypeOU(data = dat[[trait]], tree = tree, regimes = pharyngeal, 
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

p3 <- ggplot(subset(data, param %in% c("sigma"))) +
      geom_boxplot(aes(model, value, fill=regimes)) + 
      facet_wrap(trait ~ param, scales = "free_y") 

p4 <- ggplot(subset(data, param %in% c("alpha"))) +
      geom_boxplot(aes(model, value, fill=regimes)) + 
      facet_wrap(trait ~ param, scales = "free_y") 

save(list=ls(), file=sprintf("%s_lik.Rdat", id))
ggsave(sprintf("%s_lik.png", id), p1)
ggsave(sprintf("%s_params_p2.png", id),  p2)
ggsave(sprintf("%s_params_p3.png", id),  p3)
ggsave(sprintf("%s_params_p4.png", id),  p4)

#upload(sprintf("%s_*.png", id), gitaddr=gitaddr, tag="phylogenetics", save=FALSE)


