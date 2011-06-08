require(wrightscape)
require(socialR)
source("parrotfish_data.R")
# since we can't install package while other reps are running
source("../R/wrightscape.R")
source("../R/mcmc.R")
source("../R/likelihood.R")
source("../R/prior_library.R")


# from simulated_release
load("simtest.Rdat")


sfInit(parallel=T, cpu=4)
sfLibrary(wrightscape)
sfExportAll()
o <- phylo_mcmc(sim_trait, labrid$tree, intramandibular, MaxTime=1e5, indep=1e2, nchains=4)

png("alphadist.png")
cold_chain <- o$chains[[1]]
colnames(cold_chain) <- c("Pi", "Xo", "alpha1", "alpha2", "sigma", "theta")
poste_alpha1 <- density(cold_chain[, "alpha1"])
poste_alpha2 <- density(cold_chain[, "alpha2"])
xlim <- c(min(poste_alpha1$x, poste_alpha2$x), max(poste_alpha1$x, poste_alpha2$x)) 
ylim <- c(min(poste_alpha1$y, poste_alpha2$y), max(poste_alpha1$y, poste_alpha2$y)) 
plot(poste_alpha2, xlab="alpha", main="Selection Strength", xlim=xlim, ylim=ylim)
polygon(poste_alpha1, col=rgb(0,1,0,.5))
polygon(poste_alpha2, col=rgb(0,0,1,.5))
dev.off()

social_report(file="alphadist.png", tag="phylogenetics")


