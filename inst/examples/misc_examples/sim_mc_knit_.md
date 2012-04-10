# simulation.R

``` {r }
rm(list=ls())
require(wrightscape)
require(snowfall)
require(ggplot2)

# store the unique id of this script version
require(socialR)
gitaddr <- gitcommit("simulation.R")
id <- gitlog()$shortID
print(id)

data(labrids)

# rename the regimes less technically
levels(pharyngeal) = c("wrasses", "parrotfish")
regime <- pharyngeal
#  Create a dataset by simulation, where the parrotfish all have lower alpha
a1_spec  <- list(alpha = "indep", sigma = "global", theta = "global")
a1 <- multiTypeOU(data = dat["close"], tree = tree, regimes = regime, 
	     model_spec = a1_spec,  control = list(maxit=5000))

# expected "startree" standard deviations: wrasses 1.25, parrotfish 5
a1$alpha[1] <- 10
a1$alpha[2] <- 1e-10
a1$sigma <- c(5, 5)  
a1$theta <- c(0,0)   
dat[["constraint release"]] <-simulate(a1)[[1]]


# Check out the variance in the relative groups -- it should be larger in parrotfish
# in an extreme example, but need not be... 
testcase <- dat[["constraint release"]]
testcase[regime =="wrasses" & !is.na(testcase) & testcase != 0] -> lowvar
testcase[regime !="wrasses" & !is.na(testcase) & testcase != 0] -> highvar
print(c(var(lowvar), var(highvar)))


## We can repeat the whole thing with a model based on differnt sigmas, to make sure 
s1_spec  <- list(alpha = "global", sigma = "indep", theta = "global")
s1 <- multiTypeOU(data = dat["close"], tree = tree, regimes = pharyngeal, 
	     model_spec = s1_spec,  control = list(maxit=5000))
# Order of entries in sigma is the order regime names are given by levels (alphabetical)
names(s1$sigma) <- levels(regime)
s1$sigma[1] <- sqrt(2*5*1.25)
s1$sigma[2] <- sqrt(2*5*5)
s1$alpha <- c(5, 5)  # We can keep those parameters estimated from data or update them
a1$theta <- c(0,0)   
dat[["faster evolution"]] <-simulate(s1)[[1]]
testcase <- dat[["faster evolution"]]
testcase[regime == "wrasses" & !is.na(testcase) & testcase != 0 ] -> lowvar
testcase[regime != "wrasses" & !is.na(testcase) & testcase != 0 ] -> highvar
print(c(var(lowvar), var(highvar)))

## Now we have a trait where change in alpha is responsible, 
## and one in which sigma change is responsible. 
## Can we correctly identify each??

traits <- c("constraint release", "faster evolution")

sfInit(par=T, 10)    # for debugging locally
sfLibrary(wrightscape)
sfLibrary(pmc)

fits <- lapply(traits, function(trait){

	multi <- function(modelspec)
	    m <- multiTypeOU(data = dat[[trait]], tree = tree, regimes = pharyngeal, 
			     model_spec = modelspec, 
			     control = list(maxit=5000)
			    ) 

	  s1 <- multi(list(alpha = "global", sigma = "indep", theta = "global")) 
	  a1  <- multi(list(alpha = "indep", sigma = "global", theta = "global")) 
sfExportAll()

	mc <- montecarlotest(s1, a1, nboot=100)
})

````
