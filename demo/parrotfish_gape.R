# parrotfish.R
require(wrightscape)
require(pmc)
require(socialR)
tag="phylogenetics wrightscape labrids"
source("parrotfish_data.R")


general = list(alpha="indep", sigma="indep", theta="indep")
rc = list(alpha="indep", sigma="global", theta="global")


fit_input <- list(data=labrid$data["close"], tree=labrid$tree,
                  regimes=intramandibular, model_spec=rc, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50))

rc_fit <- do.call(multiTypeOU, fit_input)



save(list=ls(), file="parrotfish_gape.Rdat") #bit o checkpting
tweet("@cboettig Parrotfish run initialized")

#require(snowfall)
#sfInit(parallel=TRUE, cpu=16)
#sfLibrary(wrightscape)
#fExportAll()

#ouch_vs_gm <- montecarlotest(ouch, gm, nboot=80, cpu=16)
#social_plot(plot(ouch_vs_gm), tag=tag)


