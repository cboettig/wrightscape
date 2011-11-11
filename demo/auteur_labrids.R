
# labrid example
rm(list=ls()) # clean workspace

require(phytools)
require(auteur)
# This data has not been released
path = "../data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
diet_data <- read.csv(paste(path,"labriddata_parrotfish.csv", sep=""))

#### size correct length and weight as fraction of body mass ####
#for(i in c(3,4,6,7,8)){
#	diet_data[i] <- diet_data[i]/diet_data[5]
#}
#diet_data[5] <- log(diet_data[5]) 

corrected_data <- diet_data
# Use the simple trait abbreviations from Price et al 2010
traitnames <- c("Species", "group", "gape", "prot", "bodymass",
                "AM", "SH", "LP", "close", "open", "kt")
names(corrected_data) <- traitnames
# Lengths are log transformed 
corrected_data[["gape"]] <- log(corrected_data[["gape"]])
corrected_data[["prot"]] <- log(corrected_data[["prot"]])
#masses are log(cube-root) transformed
corrected_data[["bodymass"]] <- log(corrected_data[["bodymass"]])/3
corrected_data[["AM"]] <- log(corrected_data[["AM"]])/3
corrected_data[["SH"]] <- log(corrected_data[["SH"]])/3
corrected_data[["LP"]] <- log(corrected_data[["LP"]])/3
# ratios are fine as they are

# Drop any unmatched tip-traits
ape <- treedata(labrid_tree, corrected_data[,3:11],
                corrected_data[,1])


# Run Revell's phylogenetic size corrections
ape$data["bodysize"]
out <- phyl.resid(ape$phy, ape$data[,"bodymass"], ape$data[,c("gape", "prot","AM", "SH", "LP")] )
## phyl.resid changes order of species listing. Merge for a set of uncorrected and corrected traits.  
traits <- merge(ape$data, out$resid, by="row.names")
# columns that are transformed now have gape.x for untransformed, gape.y for transformed.  
rownames(traits) <- traits[,1]
traits <- traits[,-1]



phy <- ape$phy
dat <- traits[["prot.y"]]
names(dat) <- rownames(traits)
## run two short reversible-jump Markov chains
# (create some random strings for temporary file names)
r=paste(sample(letters,9,replace=TRUE),collapse="")

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








