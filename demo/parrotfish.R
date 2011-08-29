# parrotfish.R
rm(list=ls())

RLIBS="~/R/x86_64-redhat-linux-gnu-library/2.13"
.libPaths(c(RLIBS, .libPaths()))

require(wrightscape)
require(pmc)

############ Notebook logging header ##############
require(socialR)
script <- "parrotfish.R"
tags="phylogenetics wrightscape labrids"
gitopts <- list(user = "cboettig", dir = "demo", repo = "wrightscape") 
gitaddr <- gitcommit(script, gitopts)
####################################################

source("parrotfish_data.R")

alphas <- multiTypeOU(data=labrid$data["prot.y"], tree=labrid$tree, 
                  regimes=intramandibular, 
                  model_spec=list(alpha="indep", sigma="global", 
                  theta="global"), 
                  Xo=NULL, alpha = .1, sigma = 40, theta=NULL,
                  method ="SANN", control=list(maxit=100000,temp=50,tmax=20))

sigmas <- multiTypeOU(data=labrid$data["prot.y"], tree=labrid$tree, regimes=intramandibular, 
                model_spec=list(alpha="global", sigma="global", theta="global"), 
                  Xo=NULL, alpha = .1, sigma = 40, theta=NULL,
                  method ="SANN", control=list(maxit=100000,temp=50,tmax=20))

# models we are testing
A <- sigmas 
B <- alphas

nboot <- 64*4
require(snow)
cluster <- makeCluster(64, type="MPI")
clusterEvalQ(cluster, library(pmc))
clusterEvalQ(cluster, library(wrightscape))
clusterExport(cluster, "A")
clusterExport(cluster, "B")
A_sim <- parLapply(cluster, 1:nboot, function(x) 
                   compare_models(A, B))
B_sim <- parLapply(cluster, 1:nboot, function(x) 
                   compare_models(B, A))
stopCluster(cluster)
mc <- collect(A_sim, B_sim, A, B)
save(list=ls(), file="parrotfish.Rdat")



### commands reporting this run from head node; wll link to this code ###
# require(socialR); require(wrightscape); load("parrotfish.Rdat")
# png("parrotfish.png"); plot(mc); dev.off()
# upload("parrotfish.png", script=script, gitaddr=gitaddr, tags=tags)



