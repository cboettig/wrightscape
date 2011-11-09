# ouwie_parrotfish.R
# Author: Carl Boettiger <cboettig@gmail.com>
# Date: 19 October 2011

#==============================================================#
# Load libraries and data                                      #
#==============================================================#
rm(list=ls()) # clean workspace
require(phytools)
require(geiger)

# This data has not been released
path = "../data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
diet_data <- read.csv(paste(path,"labriddata_parrotfish.csv", sep=""))

#==============================================================#
# Size-correct the data                                        #
#==============================================================#
corrected_data <- diet_data
# Use the simple names from Price et al 2010
traitnames <- c("Species", "group", "gape", "prot", "bodymass", "AM", "SH", "LP", "close", "open", "kt")
names(corrected_data) <- traitnames
# Lengths are log transformed 
corrected_data[["gape"]] <- log(corrected_data[["gape"]])
corrected_data[["prot"]] <- log(corrected_data[["prot"]])
# masses are log(cube-root) transformed
corrected_data[["bodymass"]] <- log(corrected_data[["bodymass"]])/3
corrected_data[["AM"]] <- log(corrected_data[["AM"]])/3
corrected_data[["SH"]] <- log(corrected_data[["SH"]])/3
corrected_data[["LP"]] <- log(corrected_data[["LP"]])/3
# natural ratios are fine as they are

# We'll get just the parrotfish data by taking only the traits for parrotfish.  
# treedata() will automatically drop the tips with wrasses from the phylogeny.  
parrotfish <- corrected_data[corrected_data[,2]=="parrotfish",]
ape <- treedata(labrid_tree, parrotfish[,3:11], parrotfish[,1])
out <- phyl.resid(ape$phy, ape$data[,"bodymass"], ape$data[,c("gape", "prot","AM", "SH", "LP")] )
traits <- merge(ape$data, out$resid, by="row.names")
# columns that are transformed now have gape.x for untransformed, gape.y for transformed, etc
# put real rownames back on and drop the column of names at the beginning
rownames(traits) <- traits[,1]
traits <- traits[,-1]


phy <- ape$phy
dat <- traits[["prot.y"]]
names(dat) <- rownames(traits)
## run two short reversible-jump Markov chains
# (create some random strings for temporary file names)
r=paste(sample(letters,9,replace=TRUE),collapse="")

# run four short MCMC chains to search for a change point in brain weight
out <- lapply(1:2, 
         function(x) rjmcmc.bm(phy=phy, dat=dat,
          ngen=100000, sample.freq=10, prob.mergesplit=0.1, simplestart=TRUE,
          prop.width=1, fileBase=paste(r,x,sep=".")))

# collect directories
dirs=dir("./",pattern=paste("BM",r,sep="."))
pool.rjmcmcsamples(base.dirs=dirs, lab=r)
 
## view contents of .rda
 load(paste(paste(r,"combined.rjmcmc",sep="."), 
      paste(r,"posteriorsamples.rda",sep="."),sep="/"))
 print(head(posteriorsamples$rates))
 print(head(posteriorsamples$rate.shifts))
  
## plot Markov sampled rates
png("parrotfish.png")
shifts.plot(phy=phy, base.dir=paste(r,"combined.rjmcmc",sep="."), burnin=0.5, legend=TRUE, edge.width=4, x.lim = c(0,60))
dev.off()
# clean-up: unlink those directories
#    unlink(dir(pattern=paste(r)),recursive=TRUE)


require(socialR)
upload("parrotfish.png", script="auteur_labrids.R", tag="phylogenetics")






