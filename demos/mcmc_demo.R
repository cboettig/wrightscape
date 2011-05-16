# mcmc_demo.R
require(wrightscape)
source("parrotfish_data.R")
o <- ws_mcmc(labrid$data[1], labrid$tree, intramandibular, alpha=2, sigma=2, MaxTime=1e7)
