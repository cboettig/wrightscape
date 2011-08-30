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
hash <- gitlog()$commitID
print(hash)
####################################################

source("parrotfish_data.R")

args <- function(spec){
	  list(data=labrid$data["prot.y"], tree=labrid$tree, 
               regimes=intramandibular, model_spec = spec,                 
               Xo=NULL, alpha = .1, sigma = 1, theta=NULL,
               method ="SANN", control=list(maxit=100000,temp=50,tmax=20))
	  }

alphas <- do.call(multiTypeOU, 
          args(list(alpha="indep", sigma="global", theta="global")))

#ou <- hansen(data=labrid$data["prot.y"], tree=labrid$tree, 
#	     labrid$noregimes, sqrt.alpha=sqrt(.1), sigma=1)

# brownie hypothesis, clearly not winning (by likelihood score alone)
sigmas <- do.call(multiTypeOU, 
	   args(list(alpha="fixed", sigma="indep", theta="global")))

# probably same as alphas 
sigmas2 <- do.call(multiTypeOU, 
	   args(list(alpha="global", sigma="indep", theta="global")))

 worth attempting?
indeps <- do.call(multiTypeOU, 
          args(list(alpha="indep", sigma="indep", theta="indep")))


# models we are testing
A <- sigmas2
B <- indeps

nboot <- 64*4
require(snow)
cluster <- makeCluster(128, type="MPI")
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
save(list=ls(), file=paste(hash, ".Rdat"))



### commands reporting this run from head node; wll link to this code ###
# require(socialR); require(wrightscape); load(paste(hash, ".Rdat"))
# png("parrotfish.png"); plot(mc); dev.off()
# png("Apars.png"); plot_pars(mc$null_par_dist); dev.off()
# png("Bpars.png"); plot_pars(mc$test_par_dist); dev.off()
# upload(c("parrotfish.png", "Apars.png", "Bpars.png"),script,git=gitaddr,tag=tags)



