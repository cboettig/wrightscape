# primates.R
## load the packages we'll need
RLIBS="~/R/x86_64-redhat-linux-gnu-library/2.13"
.libPaths(c(RLIBS, .libPaths()))


require(wrightscape)
require(auteur)  # for the data
require(maticce) # to generate paintings



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
png("monkeytree.png", 2000, 2000)
plot(monkey$tree, regimes=new_world, cex=.5)
dev.off()


### Time to use wrightscape.  Note that the general model form 
# is specified by model_spec list.  This specifies which parameters
# out of alpha, theta, and sigma are independently estimated on 
# each regime, kept global across regimes, or, in the case of alpha,
# fixed to zero (to give purely Brownian behavior).  i.e.
# ouch model is equivalent to: list(alpha="global", sigma="global",
# theta="indep"), while the brownie model is equivalent to 
# list(alpha="fixed", sigma="indep", theta="global") 
# Starting parameter estimates are given or left to NULL.  
# method refers to the optimization, using parameters given in control.  
# We'll use simulated annealing to get a robust result.  


#####  Estimate the models by maximum likelihood, as in OUCH #####
bm <- brown(data=monkey$data, tree=monkey$tree)
ou <- hansen(data=monkey$data, tree=monkey$tree,
             regimes=monkey$noregimes, sigma=1, sqrt.alpha=1 )
ouch <- hansen(data=monkey$data, tree=monkey$tree,
             regimes=new_world, sigma=1, sqrt.alpha=1 )

alphas <- multiTypeOU(data=monkey$data, tree=monkey$tree,
                      regimes=new_world,model_spec = 
                      list(alpha="indep",sigma="global", 
                      theta="global"), Xo=NULL, alpha = .1, 
                      sigma = 40, theta=NULL, method ="SANN", 
                      control = list(maxit=100000,temp=50,tmax=20))

sigmas <- multiTypeOU(data=monkey$data, tree=monkey$tree, 
                      regimes=new_world, model_spec= 
                      list(alpha="fixed", sigma="indep", 
                      theta="global"), Xo=NULL, alpha = .1,
                      sigma = 40, theta=NULL, method ="SANN",
                      control=list(maxit=100000,temp=50,tmax=20))

nboot <- 64*4

require(snow)
cluster <- makeCluster(64, type="MPI")
clusterEvalQ(cluster, library(pmc))
clusterEvalQ(cluster, library(wrightscape))
clusterExport(cluster, "sigmas")
clusterExport(cluster, "alphas")
A_sim <- parLapply(cluster, 1:nboot, function(x) 
                   compare_models(sigmas,alphas))
B_sim <- parLapply(cluster, 1:nboot, function(x) 
                   compare_models(alphas, sigmas))
stopCluster(cluster)

sigmas_v_alphas <- collect(A_sim, B_sim, sigmas, alphas)
save(list=ls(), file="primates2.Rdat")




