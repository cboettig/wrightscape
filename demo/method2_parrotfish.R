# ouwie_parrotfish.R
# Author: Carl Boettiger <cboettig@gmail.com>
# Date: 19 October 2011


#=============================================================================#
# Load libraries and data                                                     #
#=============================================================================#
require(phytools)
require(geiger)

# This data has not been released
path = "../data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
diet_data <- read.csv(paste(path,"labriddata_parrotfish.csv", sep=""))


#=============================================================================#
# Size-correct the data                                                       #
#=============================================================================#
corrected_data <- diet_data
# Use the simple names from Price et al 2010
traitnames <- c("Species", "group", "gape", "prot", "bodymass",
                "AM", "SH", "LP", "close", "open", "kt")
names(corrected_data) <- traitnames
# Lengths are log transformed 
corrected_data[["gape"]] <- log(corrected_data[["gape"]])
corrected_data[["prot"]] <- log(corrected_data[["prot"]])
# masses are log(cube-root) transformed
corrected_data[["bodymass"]] <- log(corrected_data[["bodymass"]])/3
corrected_data[["AM"]] <- log(corrected_data[["AM"]])/3
corrected_data[["SH"]] <- log(corrected_data[["SH"]])/3
corrected_data[["LP"]] <- log(corrected_data[["LP"]])/3
# ratios are fine as they are

# We'll get just the parrotfish data by taking only the traits for parrotfish.  
# treedata() will automatically drop the tips with wrasses from the phylogeny.  
parrotfish <- corrected_data[corrected_data[,2]=="parrotfish",]
ape <- treedata(labrid_tree, parrotfish[,3:11], parrotfish[,1])

out <- phyl.resid(ape$phy, ape$data[,"bodymass"],
                  ape$data[,c("gape", "prot","AM", "SH", "LP")] )
## phyl.resid changes order of species listing. Merge for a set of uncorrected and corrected traits.  
traits <- merge(ape$data, out$resid, by="row.names")
# columns that are transformed now have gape.x for untransformed, gape.y for transformed.  


#================================================================================#
# Assign regimes                                                                 #
#================================================================================#


#' A function that returns the edge numbers of all branches descended from
#' the mrca of sp1 and sp2.  
#' @param phy a phylogenetic tree of class phylo (ape)
#' @param sp1 the tip.label of the first species on the tree to be used
#' @param sp1 the tip.label of the second species on the tree to be used
get_clade <- function(phy, sp1, sp2){
  M <- mrca(phy)
  who <- M[sp1, sp2]
  desc <- which(phy$edge[,1] %in% who)
  l1 <- length(desc)
  l2 <- length(desc)+1
  while(l1!=l2){
    l1 <- length(desc)
    who <- phy$edge[desc, 2]
    desc <- unique(c(desc, which(phy$edge[,1] %in% who)))
    l2 <- length(desc)
  }
  desc
}


# a list of the descendants from that common ancestor
desc <- get_clade(ape$phy, "Chlorurus_sordidus", "Hipposcarus_longiceps")
regimes <- rep(1, length(ape$phy$edge[,1])) 
regimes[desc] <- 0

# use to define regimes
colors <- rep("black", length(ape$phy$edge[,1]))
colors[desc] <- "red"
plot(ape$phy, edge.color=colors)

# pick out tips and nodes
tips <-  !(ape$phy$edge[,2] %in% ape$phy$edge[,1]) 
reg_tips <- regimes[tips]
reg_nodes <- c(1, regimes[!tips])

ape$phy$node.label <- as.character(reg_nodes)

# coloring just for fun/to check accuracy
plot(ape$phy) 
nodelabels(pch=21, bg=c("black", colors[!tips]))
tiplabels(pch=21, bg=colors[tips])

#=============================================================================#
# Data into OWie format                                                       #
#=============================================================================#

# And now for the data in OUwie format: (kinda ouch, kinda ape)
# Ape format ordering is arbitrary, horrible idea.  
# Oh well, let's match by names
parrotfish_traits <- data.frame(Reg=as.integer(reg_tips), row.names=ape$phy$tip.label)
rownames(traits) <- traits$Row.names
trait_data <- merge(parrotfish_traits, traits[,-1], by="row.names")
names(trait_data)[1] <- "Genus_species"

# Get just the active trait # not sure if necessary
parrotfish_trait <- trait_data[c("Genus_species", "Reg", "prot.y")]

paste("starting OUwie")

require(OUwie)
oumv  <- OUwie(ape$phy, parrotfish_trait, model = c("OUMV"),  
               root.station=TRUE, plot.resid=FALSE)
ouma  <- OUwie(ape$phy, parrotfish_trait, model = c("OUMA"),
               root.station=TRUE, plot.resid=FALSE)
oumva <- OUwie(ape$phy, parrotfish_trait, model = c("OUMVA"),
               root.station=TRUE, plot.resid=FALSE)






### Check that order doesn't matter (tree listing labels in different order than trait)
#data(tworegime)
#owie1 <- OUwie(tree,trait,model=c("OUMV"),root.station=FALSE, plot.resid=TRUE)
#trait <- trait[with(trait, order(X)),]
#owie2 <- OUwie(tree,trait,model=c("OUMV"),root.station=FALSE, plot.resid=TRUE)
#identical(owie1$Param.est, owie2$Param.est) # Yup, checks out
#


