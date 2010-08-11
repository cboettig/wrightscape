# labrid example

cpu=2
nboot=2

require(wrightscape)
# This data has not been released
#path = "~/data/labrids/"
path = "/home/cboettig/Documents/ucdavis/research/phylotrees/data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
fin_data <-read.csv(paste(path,"labrid.csv", sep=""))
diet_data <- read.csv(paste(path,"labriddata_parrotfish.csv", sep=""))

# size correct length and weight by body mass
for(i in c(3,4,6,7,8)){
	diet_data[i] <- diet_data[i]/diet_data[5]
}
diet_data[5] <- log(diet_data[5]) 



labrid <- format_data(labrid_tree, diet_data, species_names=diet_data[,1],  regimes = 2)  

#for(i in 3:11)
sfInit(parallel=TRUE, cpu=cpu)
sfExportAll()
sfLibrary(wrightscape)

i <- 3

#sfSapply(c(3,4,5,9,10,11), function(i){
	trait_name <- names(labrid$data)[i]	
	trait <- labrid$data[i]

	bm <- brown(trait, labrid$tree)
	ou1 <- hansen(trait, labrid$tree, regime=labrid$noregimes, .01, .01)
	ou2 <- hansen(trait, labrid$tree, regime=labrid$regimes, .01, .01)
	ws2 <- wrightscape(trait, labrid$tree, regime=labrid$regimes, (ou2@sqrt.alpha)^2, ou2@sigma, theta=ou2@theta[[1]])
#	ws2 <- wrightscape(trait, labrid$tree, regime=labrid$regimes, 1, 1, theta=ou2@theta[[1]])

	models <- list(bm = bm, ou1=ou1, ou2=ou2, ws2=ws2)

	par_boots <- fast_boot(ws2, nboot = nboot, cpus=cpu)
	png(paste(trait_name, "_pars.png", sep=""), width=800, height=400)
	plot(par_boots)
	legend("topright", levels(ws2$regimes), lty=c(1,2), lwd=2) 
	dev.off()


	LR <- choose_model(models, nboot=nboot, cpus=cpu)
		png(paste(trait_name, "_model.png", sep=""), width=2000, height=600) 
		par(mfrow=c(1,3))
		pretty_plot(LR[[1]], main="support for OU over BM")
		pretty_plot(LR[[2]], main="support for 2 peaks over 1")
		pretty_plot(LR[[3]], main="support for differential selective strength")
	dev.off()
})








