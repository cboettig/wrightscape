# parrotfish.R
require(wrightscape)
require(pmc)
require(socialR)
tag="phylogenetics wrightscape labrids"

source("parrotfish_data.R")
source("loop_models_traits_regimes.R")

model_list <- list("brown", "hansen", "ouch", "brownie", "wright", "release_constraint")
regime_list <-  list(intramandibular=intramandibular)

test <- fit_all(model_list[c(4,6)], labrid$data[10:13], regime_list, labrid$tree)

llik_matrix(test, 1)
alpha_matrix(test, trait=1)




## sometimes hansen returns a silly large value of alpha, which will break the generalized likelihood function (too stiff)
 #h <- fit("hansen", labrid$data[13], intramandibular, labrid$tree, .01, .01)
 # w <- fit("release_constraint", labrid$data[13], intramandibular, labrid$tree, h@sqrt.alpha^2, h@sigma)


