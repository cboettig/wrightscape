# parrotfish.R
rm(list=ls())
require(wrightscape)
require(warningsignals)



############ Notebook logging header ##############
require(socialR)
script <- "fast_ML_ex.R"
tags="phylogenetics wrightscape labrids"
gitopts <- list(user = "cboettig", dir = "demo", repo = "wrightscape") 
gitaddr <- gitcommit(script, gitopts)
on.exit(system("git push")) 
tags <- "phylogenetics"  
tweet_errors(script, gitopts, tags) 
####################################################

source("parrotfish_data.R")

alphas <- multiTypeOU(data=labrid$data["close"], tree=labrid$tree, 
                  regimes=intramandibular, 
                  model_spec=list(alpha="indep", sigma="global", 
                  theta="global"), 
                  Xo=NULL, alpha = .1, sigma = 10, theta=NULL)

sigmas <- multiTypeOU(data=labrid$data["close"], tree=labrid$tree, regimes=intramandibular, 
                model_spec=list(alpha="fixed", sigma="indep", theta="global"), 
                  Xo=NULL, alpha = .1, sigma = 10, theta=NULL)

bm <- brown(data=labrid$data["close"], tree=labrid$tree)
ou <- hansen(data=labrid$data["close"], tree=labrid$tree, regimes=intramandibular, 
             sqrt.alpha = .1, sigma = 10)



require(snowfall)
sfInit(parallel=TRUE, cpu=16)
sfLibrary(wrightscape)
sfLibrary(warningsignals)
sfLibrary(ouch)
sfExportAll()

boots <- montecarlotest(bm, ou, nboot=20, cpu=16)
#boots <- montecarlotest(sigmas, alphas, nboot=200, cpu=8)
png("sigmas_v_alphas.png")
plot(boots)
dev.off()
upload("sigmas_v_alphas.png", script=script, tags=tags, gitaddr=gitaddr)


