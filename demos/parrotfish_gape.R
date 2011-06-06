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

test <- wright(data=labrid$data["gape.y"], tree=labrid$tree,
                  regimes=intramandibular, alpha = .1, sigma = .1)


rc <- do.call(release_constraint,fit_input)
gm <- do.call(wright, fit_input)
multi_peak <- do.call(hansen, fit_input)
bm <- do.call(bm, fit_input)
brownie <- do.call(brownie, fit_input)


save(list=ls(), file="parrotfish_gape.Rdat") #bit o checkpting
brownie_vs_gm <- montecarlotest(brownie, gm, nboot=80, cpu=16)
save(list=ls(), file="parrotfish_gape.Rdat")
multi_peak_vs_gm montecarlotest(mulit_peak, gm, nboot=80, cpu=16)
save(list=ls(), file="parrotfish_gape.Rdat")



