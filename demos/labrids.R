# labrid example

require(wrightscape)
# This data has not been released
path = "/home/cboettig/Documents/ucdavis/research/phylotrees/data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
labrid_data <-read.csv(paste(path,"labrid.csv", sep=""))
diet_data <- read.csv(paste(path,"labrid_herbivory.csv", sep=""))

labrid <- format_data(labrid_tree, diet_data, species_names=diet_data[,1],  regimes = 2)  

gapesize <- labrid$data[3]/labrid$data[5]
bm <- brown(gapesize, labrid$tree)

ou1 <- hansen(gapesize, labrid$tree, regime=labrid$noregimes, 1, 1)
ws1 <- wrightscape(gapesize, labrid$tree, regime=labrid$noregimes, 1, 1)
ou2 <- hansen(gapesize, labrid$tree, regime=labrid$regimes, 1, 1)
ws2 <- wrightscape(gapesize, labrid$tree, regime=labrid$regimes, (ou2@sqrt.alpha)^2, ou2@sigma)


par_boots <- fast_boot(ws2, nboot = 100, cpus=16)
png("pars.png", width=800, height=400)
plot(par_boots)
legend("topright", c("Wrasses", "Parrotfish"), lty=c(1,2), lwd=2) 
dev.off()


labrid_models <- list(bm = bm, ws1=ws1, ou2=ou2, ws2=ws2)
LR <- choose_model(labrid_models, 100, cpus=16)
	png("model.png",width=2000, height=600) 
	par(mfrow=c(1,3))
	pretty_plot(LR[[1]], main="support for OU over BM")
	pretty_plot(LR[[2]], main="support for 2 peaks over 1")
	pretty_plot(LR[[3]], main="support for differential selective strength")
dev.off()









