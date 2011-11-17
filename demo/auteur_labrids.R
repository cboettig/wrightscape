
# labrid example
rm(list=ls()) # clean workspace

require(auteur)
require(wrightscape)
data(labrids)

phy <- ape.phy
dat <- ape.traits[["prot.y"]]
names(dat) <- rownames(traits)


## run two short reversible-jump Markov chains
# (create some random strings for temporary file names)
r=paste(sample(letters,9,replace=TRUE),collapse="")

# run four short MCMC chains to search for a change point in brain weight
## Get ready for some MPI parallel computing action
require(snow)
cl <- makeCluster(20, type="MPI")
clusterExport(cl, c("phy", "dat", "r"))
clusterEvalQ(cl, library(auteur))
## RUN AUTEUR
out <- parLapply(cl, 1:20, 
         function(x) rjmcmc.bm(phy=phy, dat=dat,
          ngen=1000000, sample.freq=10, prob.mergesplit=0.1, simplestart=TRUE,
          prop.width=1, fileBase=paste(r,x,sep=".")))
stopCluster(cl)

# collect directories
dirs=dir("./",pattern=paste("BM",r,sep="."))
pool.rjmcmcsamples(base.dirs=dirs, lab=r)
 
## view contents of .rda
 load(paste(paste(r,"combined.rjmcmc",sep="."), 
      paste(r,"posteriorsamples.rda",sep="."),sep="/"))
 print(head(posteriorsamples$rates))
 print(head(posteriorsamples$rate.shifts))
 
save(list=ls(), file="auteur_labrids.Rdat")
 
## plot Markov sampled rates
cairo_pdf("auteur_labrids.pdf")
shifts.plot(phy=phy, base.dir=paste(r,"combined.rjmcmc",sep="."), burnin=0.5, legend=TRUE, edge.width=4, x.lim = c(0,60))
dev.off()



#require(socialR)
#upload("labrids.png", script="auteur_labrids.R", tag="phylogenetics")

# clean-up: unlink those directories
#    unlink(dir(pattern=paste(r)),recursive=TRUE)








