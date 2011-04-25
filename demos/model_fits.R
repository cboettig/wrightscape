# labrid example

require(wrightscape)
require(pmc)
require(socialR)

# This data has not been released
path = "../../../data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
fin_data <-read.csv(paste(path,"labrid.csv", sep=""))
diet_data <- read.csv(paste(path,"labriddata_parrotfish.csv", sep=""))

# size correct length and weight by body mass
for(i in c(3,4,6,7,8)){
	diet_data[i] <- diet_data[i]/diet_data[5]
}
diet_data[5] <- log(diet_data[5]) 
labrid <- format_data(labrid_tree, diet_data, species_names=diet_data[,1],  regimes = 2)  

#for(i in 3:11){


#sfInit(parallel=TRUE, cpu=2)
#sfExportAll()
#sfLibrary(wrightscape)
#sfLibrary(pmc)
i <- 3

#sfSapply(c(3,4,5,9,10,11), function(i){
	trait_name <- names(labrid$data)[i]	
	trait <- labrid$data[i]

	bm <- brown(trait, labrid$tree)
	ou1 <- hansen(trait, labrid$tree, regime=labrid$noregimes, .01, .01)
	ou2 <- hansen(trait, labrid$tree, regime=labrid$regimes, .01, .01)
	ws2 <- wrightscape(trait, labrid$tree, regime=labrid$regimes, (ou2@sqrt.alpha)^2, ou2@sigma, theta=ou2@theta[[1]])

    ouch_test <- ouch(trait, labrid$tree, regime=labrid$regimes)
    brownie_test <- brownie(trait, labrid$tree, regime=labrid$regimes)

#})



