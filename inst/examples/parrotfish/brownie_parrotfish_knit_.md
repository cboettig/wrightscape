# ouwie_parrotfish.R
# Author: Carl Boettiger <cboettig@gmail.com>
# Date: 19 October 2011

#==============================================================#
# Load libraries and data                                      #
#==============================================================#
rm(list=ls()) # clean workspace
require(phytools)
require(geiger)
require(OUwie)

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



#==============================================================#
# Paint regimes onto the phylogeny                             #
#==============================================================#
source("method2_tools.R")

clade <- get_clade(ape$phy, "Chlorurus_sordidus", "Hipposcarus_longiceps")
tree <- ape$phy
tree$mapped.edge<-matrix(0,nrow(tree$edge),2, dimnames=list(edge=apply(tree$edge,1,function(x) paste(x,collapse=",")),state=c("0","1")))
w <- tree$edge[,2] %in% clade
tree$mapped.edge[w, "1"] <- tree$edge.length[w]  
tree$mapped.edge[!w, "0"] <- tree$edge.length[!w]  

bm <- fitContinuous(tree, traits["open"])
bml <- brownie.lite(tree, traits["open"]) 
#input <- paint_phy(ape$phy, traits,  c("Chlorurus_sordidus", "Hipposcarus_longiceps"))







save(list=ls(), file="brownie_parrotfish.Rdat")


