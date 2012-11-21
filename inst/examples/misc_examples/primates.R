# File: primates.R
# Author: Carl Boettiger <cboettig@gmail.com>
# License: BSD
# Date: 2011-11-16

rm(list=ls()) # start clean

## load the packages we'll need
require(wrightscape)
require(auteur)  # for the data

## Load and format the data 
data(primates)
monkey <- format_data(primates$phy, primates$dat)

## Paint the tree with a transition at New World monkeys
new_world_ancestor <- mrcaOUCH(c("Ateles_belzebuth",
                                 "Leontopithecus_caissara"), 
                               monkey$tree)
new_world <- paintBranches(new_world_ancestor, monkey$tree,
                            c("OldWorld", "NewWorld"))

## take a quick look at the tree
#plot(monkey$tree, regimes=new_world, cex=.5)

#####  Estimate the models by maximum likelihood, as in OUCH #####
alphas <- multiTypeOU(data=monkey$data, tree=monkey$tree,
                      regimes=new_world,model_spec = 
                      list(alpha="indep",sigma="global", 
                      theta="indep"), control=list(maxit=2000))

sigmas <- multiTypeOU(data=monkey$data, tree=monkey$tree,
                      regimes=new_world,model_spec = 
                      list(alpha="fixed",sigma="indep", 
                      theta="global"))

full <- multiTypeOU(data=monkey$data, tree=monkey$tree,
                      regimes=new_world, control=list(maxit=5000))
#                      method ="SANN", 
#                      control = list(maxit=100000,temp=50,tmax=20))


## snow cluster run
nboot <- 64
require(snow)
cluster <- makeCluster(64, type="MPI")  # for supercomputers
#cluster <- makeCluster(2, type="SOCK")    # for debugging locally

clusterEvalQ(cluster, library(wrightscape))
clusterExport(cluster, "full") # can just export ls()
full_boot <- parSapply(cluster, 1:nboot, function(x) bootstrap(full))
stopCluster(cluster)

summary(full, full_boot)
save(list=ls(), file="primates.Rdat")



## simple bootstrap
#require(snowfall); sfInit(par=T, cpu=4); sfExport("alphas"); sfLibrary(wrightscape)
#alphas_boot <- sfSapply(1:4, function(i) bootstrap(alphas))
#summary(alphas, alphas_boot)

