# intro_example.R
require(wrightscape)
require(geiger)

######## STEP 1: Load all the data ###################
path = "../data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
diet_data <- read.csv(paste(path,"labriddata_parrotfish.csv", sep=""))
# Make a copy of the data that we'll apply transformations to.  
corrected_data <- diet_data
# Use the simple names from Price et al 2010
traitnames <- c("Species", "group", "gape", "prot", "bodymass", "AM", "SH", "LP", "close", "open", "kt")
names(corrected_data) <- traitnames



###### Step 2: Transform data appropriately ######################

# Lengths are log transformed 
corrected_data[["gape"]] <- log(corrected_data[["gape"]])
corrected_data[["prot"]] <- log(corrected_data[["prot"]])
#masses are log(cube-root) transformed
corrected_data[["bodymass"]] <- log(corrected_data[["bodymass"]])/3
corrected_data[["AM"]] <- log(corrected_data[["AM"]])/3
corrected_data[["SH"]] <- log(corrected_data[["SH"]])/3
corrected_data[["LP"]] <- log(corrected_data[["LP"]])/3
# ratios are fine as they are

# We'll get just the parrotfish data by taking only the traits for parrotfish.  
# treedata() will automatically drop the tips with wrasses from the phylogeny.  
parrotfish <- corrected_data[corrected_data[,2]=="parrotfish",]
ape <- treedata(labrid_tree, parrotfish[,3:11], parrotfish[,1])


## Use Liam's package to do SIZE CORRECTION 
require(phytools)
ape$data["bodysize"]
out <- phyl.resid(ape$phy, ape$data[,"bodymass"], ape$data[,c("gape", "prot","AM", "SH", "LP")] )
## phyl.resid changes order of species listing. Merge for a set of uncorrected and corrected traits. #annoying
traits <- merge(ape$data, out$resid, by="row.names")
# columns that are transformed now have gape.x for untransformed, gape.y for transformed, etc.....  

# format data gets regimes from column specified in "regimes" (e.g. this is a column id, not a # of regimes)
# This also converts the tree and data into ouch format
labrid <- format_data(labrid_tree, traits[,2:length(traits)], species_names=traits[,1])  



################## Step 3:   PAINTING REGIMES #############
## Having taken care of traits, we paint on the 3 regime models.  
# Select common ancestor of a Chlorurus and a Hipposcarus as the changepoint
intra_ancestor <- mrcaOUCH(c("Chlorurus_sordidus", "Hipposcarus_longiceps"), labrid$tree)
intramandibular <- paintBranches(intra_ancestor, labrid$tree, c("other","intramandibular"))




##### Step 4: estimate the models, as in OUCH ############

alphas <- multiTypeOU(data=labrid$data["close"], tree=labrid$tree, 
                  regimes=intramandibular, 
                  model_spec=list(alpha="indep", sigma="global", 
                  theta="global"), 
                  Xo=NULL, alpha = .1, sigma = 40, theta=NULL,
                  method ="SANN", control=list(maxit=100000,temp=50,tmax=20))

sigmas <- multiTypeOU(data=labrid$data["close"], tree=labrid$tree, regimes=intramandibular, 
                model_spec=list(alpha="fixed", sigma="indep", theta="global"), 
                  Xo=NULL, alpha = .1, sigma = 40, theta=NULL,
                  method ="SANN", control=list(maxit=100000,temp=50,tmax=20))



