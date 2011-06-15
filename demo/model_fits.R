# labrid example
require(wrightscape)

# This data has not been released
path = "../data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
fin_data <-read.csv(paste(path,"labrid.csv", sep=""))
diet_data <- read.csv(paste(path,"labriddata_parrotfish.csv", sep=""))
# size correct length and weight by body mass
for(i in c(3,4,6,7,8)){
	diet_data[i] <- diet_data[i]/diet_data[5]
}
diet_data[5] <- log(diet_data[5]) 
labrid <- format_data(labrid_tree, diet_data, species_names=diet_data[,1],  regimes = 2)  
# Select common ancestor of a Chlorurus and a Hipposcarus as the changepoint
intra_ancestor <- mrcaOUCH(c("Chlorurus_sordidus", "Hipposcarus_longiceps"), labrid$tree)
intramandibular <- paintBranches(intra_ancestor, labrid$tree, c("other","intramandibular"))
# Select common ancestor for all parrot fish
pharyngeal_ancestor<-mrcaOUCH(c("Bolbometopon_muricatum", "Sparisoma_radians"), labrid$tree)
pharyngeal <- paintBranches(pharyngeal_ancestor, labrid$tree, c("other","pharyngeal"))
two_shifts <- paintBranches(c(pharyngeal_ancestor, intra_ancestor), labrid$tree, c("wrasses", "pharyngeal", "intramandibular") )

## This leaves the branch on which the second transition occurs unspecified (fourth regime).  
## We have to fix this manually
two_shifts[as.numeric(intra_ancestor)] <- "intramandibular"
two_shifts <- as.factor(as.character(two_shifts))
names(two_shifts) <- names(intramandibular)

### We now have access to the following configurations:
## intramandibular, pharyngeal and two_shifts paintings, (and labrid$noregimes),
## and all traits and size-corrected traits in labrid$data
## and tree in labrid$tree


  ## Note that this will converge poorly with the .01, .01 starting conditions
  #  ou3 <- hansen(trait, labrid$tree, regime=two_shifts, .01, .01 )
  #  loglik(ou3) - loglik(ou2_phar)



