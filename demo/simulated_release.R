require(wrightscape)
source("parrotfish_data.R")
general = list(alpha="indep", sigma="indep", theta="indep")
rc = list(alpha="indep", sigma="global", theta="global")
fit_input <- list(data=labrid$data["close"], tree=labrid$tree,
                  regimes=intramandibular, model_spec=rc, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50))
true <- do.call(multiTypeOU, fit_input)
true$Xo <- 0
true$alpha <- c(.01, 20)
true$sigma <- 10
true$theta <- 0
sim_trait <- simulate(true)
## Now we're up and running with fake data

fit <- multiTypeOU(data=sim_trait, tree=labrid$tree,
                  regimes=intramandibular, model_spec=rc, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50))

gen_fit <- multiTypeOU(data=sim_trait, tree=labrid$tree,
                  regimes=intramandibular, model_spec=general, 
                  Xo=NULL, alpha = .1, sigma = .1, theta=NULL,
                  method ="SANN", control=list(maxit=80000,temp=25,tmax=50))

boots <- montecarlotest(fit, gen_fit, nboot=80, cpu=16)
social_plot(plot(boots), tags="phylogenetics")

sfStop()
sfInit(parallel=T, cpu=4)
sfLibrary(wrightscape)
sfExportAll()
o <- phylo_mcmc(sim_trait, labrid$tree, intramandibular, MaxTime=1e5, indep=1e2)

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



