# mcmc_demo.R
require(wrightscape)
require(socialR)
source("parrotfish_data.R")
#sfInit(parallel=T, cpu=2)
#sfLibrary(wrightscape)
#sfExportAll()

o <- general_mcmc(labrid$data[1], labrid$tree, intramandibular, alpha=2, sigma=2, MaxTime=1e3)
colnames(o) <- c("Pi", "Xo", "alpha1", "alpha2", "sigma", "theta")
require(coda)
M <- mcmc(o)
plot(M)

