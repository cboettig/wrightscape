# labrid example

require(wrightscape)
require(pmc)
require(socialR)

# This data has not been released
path = "../data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
fin_data <-read.csv(paste(path,"labrid.csv", sep=""))
diet_data <- read.csv(paste(path,"labriddata_parrotfish.csv", sep=""))

# Select common ancestor of a Chlorurus and a Hipposcarus as the changepoint
intramandibular <- paintBranches(mrcaOUCH(c("Chlorurus_sordidus", "Hipposcarus_longiceps"), labrid$tree), labrid$tree, c("other","intramandibular"))



# size correct length and weight by body mass
for(i in c(3,4,6,7,8)){
	diet_data[i] <- diet_data[i]/diet_data[5]
}
diet_data[5] <- log(diet_data[5]) 
labrid <- format_data(labrid_tree, diet_data, species_names=diet_data[,1],  regimes = 2)  

#for(i in 3:11){
i <- 3

#sfSapply(c(3,4,5,9,10,11), function(i){
	trait_name <- names(labrid$data)[i]	
	trait <- labrid$data[i]

	bm <- brown(trait, labrid$tree)
	ou1 <- hansen(trait, labrid$tree, regime=labrid$noregimes, .01, .01)
	ou2 <- hansen(trait, labrid$tree, regime=labrid$regimes, .01, .01)


## needs to be fixed, singular lik problems?
  ouch_test <- ouch(trait, labrid$tree, regime=labrid$regimes, alpha=(ou2@sqrt.alpha)^2, sigma=ou2@sigma)
	a <- simulate(ouch_test)
	# update(ouch_test, a)



  ws2_phar <- wrightscape(trait, labrid$tree, regime=labrid$regimes, (ou2@sqrt.alpha)^2, ou2@sigma, theta=ou2@theta[[1]])
	ou2 <- hansen(trait, labrid$tree, regime=intramandibular, .01, .01)
  ws2_intra <- wrightscape(trait, labrid$tree, regime=intramandibular, (ou2@sqrt.alpha)^2, ou2@sigma, theta=ou2@theta[[1]])

loglik(ws2_par)
loglik(ws2_intra)

#  pharyngeal vs intramnadibular regimes


	brownie_test <- brownie(trait, labrid$tree, regime=labrid$regimes, sigma=ou2@sigma)
	a <- simulate(brownie_test)
	b <- update(brownie_test, a)
	a <- simulate(b)

	wright_test <- wright(trait, labrid$tree, regime=labrid$regimes, alpha=(ou2@sqrt.alpha)^2, sigma=ou2@sigma)
	a <- simulate(wright_test)
	w <- update(wright_test, a)
	a <- simulate(w)

cpu <- 1
#sfInit(parallel=TRUE, cpu=cpu)
#sfExportAll()
#sfLibrary(wrightscape)
#sfLibrary(pmc)

#out <- montecarlotest(brownie_test, ws2, cpu=1,nboot=2) 
out <- montecarlotest(brownie_test, wright_test, cpu=cpu,nboot=2) 
social_plot(plot(out), tag="phylogenetics wrightscape labrids")
#})



