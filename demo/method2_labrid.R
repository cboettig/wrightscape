# labrid example
rm(list=ls()) # clean workspace

require(phytools)
require(geiger)
require(OUwie)
# This data has not been released
path = "../data/"
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
         list(c("Bolbometopon_muricatum", "Sparisoma_radians"), 
              c("Chlorurus_sordidus", "Hipposcarus_longiceps")))

# Get just the active trait # not sure if necessary
X <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")

trait <- input$data[c("Genus_species", "Reg", "bodymass")]

require(geiger)
phy <- lambdaTree(input$phy, .0001)
trait[,3] <- rnorm(length(trait[,3]), 0, sd=2)
ouma <- OUwie(input$phy, trait, model = c("OUMA"),
               root.station=TRUE, plot.resid=FALSE)

# expected variance
ouwie_var <- .5*ouma$Param.est["sigma.sq",]/ouma$Param.est["alpha",]*(1-exp(-2*ouma$Param.est["alpha",]))
print(ouwie_var)

require(wrightscape)
data(labrids)
test <- dat[["prot.y"]]
test[!is.na(test)] <- trait[,3]  
tree <- convert(lambdaTree(convert(tree), .0001))
ws <- multiTypeOU(test, tree, two_shifts, control=list(maxit=10000), model=list(alpha="indep", sigma="global", theta="indep"))

## expected variance 
ws_var <- (1-exp(-2*ws$alpha))*.5*ws$sigma^2/ws$alpha

print(ws_var)

var(test, na.rm=T)
