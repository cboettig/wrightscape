# parrotfish.R
rm(list=ls())

RLIBS="~/R/x86_64-redhat-linux-gnu-library/2.13"
.libPaths(c(RLIBS, .libPaths()))

require(wrightscape)
require(pmc)

############ Notebook logging header ##############
require(socialR)
script <- "parrotfish_lik.R"
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

sigmas <- do.call(multiTypeOU, 
	   args(list(alpha="fixed", sigma="indep", theta="global")))

sigmas2 <- do.call(multiTypeOU, 
	   args(list(alpha="global", sigma="indep", theta="global")))

indeps <- do.call(multiTypeOU, 
          args(list(alpha="indep", sigma="indep", theta="indep")))

print(alphas$loglik)

print(c(sigmas$loglik, sigmas2$loglik, indeps$loglik))


## Compare these models from wrightscape and ouch fits:
ws_ouch <- do.call(multiTypeOU, 
	   args(list(alpha="global", sigma="global", theta="indep")))

ws_ou <- do.call(multiTypeOU, 
	   args(list(alpha="global", sigma="global", theta="global")))

ws_bm <- do.call(multiTypeOU, 
	   args(list(alpha="fixed", sigma="global", theta="global")))


bm <- brown(data=labrid$data["prot.y"], tree=labrid$tree)
ouch <- hansen(data=labrid$data["prot.y"], tree=labrid$tree, 
	     regimes=intramandibular, sqrt.alpha=sqrt(.1), sigma=1)
ou <- hansen(data=labrid$data["prot.y"], tree=labrid$tree, 
	     labrid$noregimes, sqrt.alpha=sqrt(.1), sigma=1)


print(c(bm@loglik, ws_bm$loglik))
print(c(ou@loglik, ws_ou$loglik))
print(c(ouch@loglik, ws_ouch$loglik))


save(list=ls(), file="parrotfish_loglik.rda")
