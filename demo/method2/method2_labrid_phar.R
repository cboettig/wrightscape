# labrid example
rm(list=ls()) # clean workspace
require(phytools)
require(geiger)
require(OUwie)
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


#==========================#
# Data into OWie format    #
#==========================#

# a list of the descendants from that common ancestor
source("method2_tools.R")

input <- paint_phy(ape$phy, traits,  
         c("Bolbometopon_muricatum", "Sparisoma_radians"))

# Get just the active trait # not sure if necessary
X <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")

require(snowfall)
sfInit(parallel=T, cpu=16)
sfLibrary(OUwie)
sfExportAll()

all_oumva <- sfLapply(X, function(x){
  trait <- input$data[c("Genus_species", "Reg", x)]
  oumva <- OUwie(input$phy, trait, model = c("OUMVA"),
               root.station=TRUE, plot.resid=FALSE)
})
all_ouma <- sfLapply(X, function(x){
  trait <- input$data[c("Genus_species", "Reg", x)]
  ouma <- OUwie(input$phy, trait, model = c("OUMA"),
               root.station=TRUE, plot.resid=FALSE)
})
all_oumv <- sfLapply(X, function(x){
  trait <- input$data[c("Genus_species", "Reg", x)]
  oumv <- OUwie(input$phy, trait, model = c("OUMV"),
               root.station=TRUE, plot.resid=FALSE)
})


save(list=ls(), file="method2_labrid_phar.Rdat")


