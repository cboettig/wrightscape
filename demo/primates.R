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


## Pain the tree with a transition at New World monkeys
new_world_ancestor <- mrcaOUCH(c("Ateles_belzebuth",
                                 "Leontopithecus_caissara"), 
                               monkey$tree)
new_world <- paintBranches(new_world_ancestor, monkey$tree,
                            c("OldWorld", "NewWorld"))

## take a quick look at the tree
plot(monkey$tree, regimes=new_world)



#####  Estimate the models by maximum likelihood, as in OUCH #####
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


chains <- function(){ 
          phylo_mcmc(monkey$data, monkey$tree, new_world, 
                     MaxTime=1e5, alpha=.1, sigma=.1, 
                     theta=NULL, Xo=NULL, model_spec=
                     list(alpha="indep", sigma="global",theta="global"),
                     stepsizes=0.05)[[1]]
}

N <- 10 # n-1 from the threads allocated(?) 

require(Rmpi)
mpi.spawn.Rslaves(nslaves=N)

## Clean-up 
.Last <- function(){
    if (is.loaded("mpi_initialize")){
        if (mpi.comm.size(1) > 0){
            print("Please use mpi.close.Rslaves() to close slaves.")
            mpi.close.Rslaves()
        }
        print("Please use mpi.quit() to quit R")
        .Call("mpi_finalize")
    }
}

mpi.bcast.Robj2slave(monkey)
mpi.bcast.Robj2slave(new_world)
mpi.bcast.Robj2slave(chains)

mcmc_out <- mpi.remote.exec(chains())


# close slaves and exit
mpi.close.Rslaves()
mpi.quit(save = "no")


