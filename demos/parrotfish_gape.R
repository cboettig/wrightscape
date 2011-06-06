# parrotfish.R
require(wrightscape)
require(pmc)
require(socialR)
tag="phylogenetics wrightscape labrids"
source("parrotfish_data.R")
source("loop_models_traits_regimes.R")

model_list <- list("brown", "hansen", "ouch", "brownie", "wright", "release_constraint")
regime_list <-  list(intramandibular=intramandibular)

fit_input <- list(data=labrid$data["gape.y"], tree=labrid$tree,
                  regimes=intramandibular, alpha = .1, sigma = .1,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50))
brownie_input <- list(data=labrid$data["gape.y"], tree=labrid$tree,
                  regimes=intramandibular, sigma = .1,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50))

# superfast, but not nec max w/o SANN
test <- wright(data=labrid$data["gape.y"], tree=labrid$tree,
                  regimes=intramandibular, alpha = .1, sigma = .1)

## fit with SANN options above
brownie <- do.call(brownie, brownie_input)
gm <- do.call(wright, fit_input)


save(list=ls(), file="parrotfish_gape.Rdat") #bit o checkpting
tweet("@cboettig Parrotfish run initialized")

require(snowfall)
sfInit(parallel=TRUE, cpu=16)
sfLibrary(wrightscape)
sfExportAll()

brownie_vs_gm <- montecarlotest(brownie, gm, nboot=80, cpu=16)
save(list=ls(), file="parrotfish_gape.Rdat")
social_plot(plot(brownie_vs_gm), tag=tag)

#multi_peak_vs_gm montecarlotest(mulit_peak, gm, nboot=80, cpu=16)
#save(list=ls(), file="parrotfish_gape.Rdat")
#social_plot(plot(multi_peak_vs_gm), tag=tag)


