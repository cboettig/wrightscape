# mcmc_demo.R
require(wrightscape)
require(socialR)
tags <- c("phylogenetics mcmcmc parrotfish")
source("parrotfish_data.R")
sfInit(parallel=T, cpu=4)
sfLibrary(wrightscape)
sfExportAll()

o <- general_mcmc(labrid$data[1], labrid$tree, intramandibular,
                  alpha=2, sigma=2, MaxTime=1e5, indep=1e1, stepsizes=.2)

png("convergenceTemp.png")
  burnin <- 1:1e4
  plot(o[[1]][-burnin,1], type="l")
  lines(o[[2]][-burnin,1], col="blue")
  lines(o[[3]][-burnin,1], col="green")
  lines(o[[4]][-burnin,1], col="red")
dev.off()

social_report(file="convergenceTemp.png", tag=tags, comment="MaxTime=1e5, indep=1e1, stepsizes=.2")

social_plot(hist(o[[1]][,2]), tag=tags)

require(coda)
#M <- mcmc(o[[1]])
#plot(M)

