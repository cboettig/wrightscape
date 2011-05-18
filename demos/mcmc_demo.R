# mcmc_demo.R
require(wrightscape)
require(socialR)
tags <- c("phylogenetics mcmcmc parrotfish")
source("parrotfish_data.R")
#sfInit(parallel=T, cpu=2)
#sfLibrary(wrightscape)
#sfExportAll()

o <- general_mcmc(labrid$data[1], labrid$tree, intramandibular,
                  alpha=2, sigma=2, MaxTime=1e6, indep=1e3)
colnames(o[[1]]) <- c("Pi", "Xo", "alpha1", "alpha2", "sigma", "theta")

png("convergenceTemp.png")
  plot(o[[1]][,1], type="l")
  lines(o[[2]][,1], col="blue")
  lines(o[[3]][,1], col="green")
  lines(o[[4]][,1], col="red")
dev.off()

social_report(file="convergenceTemp.png", tag=tags)

social_plot(hist(o[[1]][,2]), tag=tags)

require(coda)
#M <- mcmc(o[[1]])
#plot(M)

