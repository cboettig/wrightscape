# labrid example

require(wrightscape)
require(pmc)
require(socialR)

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
intramandibular <- paintBranches(mrcaOUCH(c("Chlorurus_sordidus", "Hipposcarus_longiceps"), labrid$tree), labrid$tree, c("other","intramandibular"))



#for(i in 3:11){
i <- 3

#sfSapply(c(3,4,5,9,10,11), function(i){
	trait_name <- names(labrid$data)[i]	
	trait <- labrid$data[i]

	bm <- brown(trait, labrid$tree)
	ou1 <- hansen(trait, labrid$tree, regime=labrid$noregimes, .01, .01)
	ou2_phar <- hansen(trait, labrid$tree, regime=labrid$regimes, .01, .01)
	ou2_intra <- hansen(trait, labrid$tree, regime=intramandibular, .01, .01)


## needs to be fixed, singular lik problems?
  ouch_phar <- ouch(trait, labrid$tree, regime=labrid$regimes, alpha=(ou2_phar@sqrt.alpha)^2, sigma=ou2_phar@sigma)
	a <- simulate(ouch_phar)
	# update(ouch_phar, a)


  ## Compare intramandibular and pharyngeal joint paintings
  ws2_phar <- wrightscape(trait, labrid$tree, regime=labrid$regimes, (ou2_phar@sqrt.alpha)^2, ou2_phar@sigma, theta=ou2_phar@theta[[1]])
  ws2_intra <- wrightscape(trait, labrid$tree, regime=intramandibular, (ou2_intra@sqrt.alpha)^2, ou2_intra@sigma, theta=ou2_intra@theta[[1]])

  loglik(ws2_phar)
  loglik(ws2_intra)

#  pharyngeal vs intramnadibular regimes
	brownie_test <- brownie(trait, labrid$tree, regime=labrid$regimes, sigma=ou2_phar@sigma)
	wright_test <- wright(trait, labrid$tree, regime=labrid$regimes, alpha=(ou2_phar@sqrt.alpha)^2, sigma=ou2_phar@sigma)

	brownie_intra <- brownie(trait, labrid$tree, regime=intramandibular, sigma=ou2_intra@sigma)
	wright_intra <- wright(trait, labrid$tree, regime=intramandibular, alpha=(ou2_intra@sqrt.alpha)^2, sigma=ou2_intra@sigma)






cpu <- 16
nboot <- 160
sfInit(parallel=TRUE, cpu=cpu)
sfExportAll()
sfLibrary(wrightscape)  # need all this just to export wrightscape?
sfLibrary(pmc)

out <- montecarlotest(brownie_test, ws2, cpu=cpu,nboot=nboot) 
social_plot(plot(out), tag="phylogenetics wrightscape labrids ws2")
out2 <- montecarlotest(brownie_test, wright_test, cpu=cpu,nboot=nboot) 
social_plot(plot(out2), tag="phylogenetics wrightscape labrids")
#})



