# carnivora.2R
nboot <- 100
cpu = 8

require(wrightscape)
path = "~/data/carnivora/"


# uses marine, freshwater, terrestrial as regimes
tree <- read.nexus(paste(path, "carnivora.nex", sep=""))
data <- read.csv(paste(path, "carnivora.txt", sep=""), sep="\t")
out <- format_data(tree, data, species_names = data[,1], regimes = data[,4])

size <- log(out$data[,2])
names(size) <- rownames(out$data)

bm <- brown(size, out$tree )
ou1 <- hansen(size, out$tree, regimes=out$noregimes, 1,1 )
ou2 <- hansen(size, out$tree, regimes=out$regimes, 1,1 )
ws2 <- wrightscape(size, out$tree, regime=out$regimes, 1,1)

models <- list(bm = bm, ou1 = ou1, ou2 = ou2, ws2 = ws2)
LR <- choose_model(models, nboot=nboot, cpu=cpu)
par_boots <- fast_boot(ws2, nboot = nboot, cpu=cpu)

png("carnivora_models_mft.png",width=1500, height=600) 
	par(mfrow=c(1,(length(models)-1) ))
	for(i in 1:(length(models)-1) ){
		pretty_plot(LR[[i]], main=paste("support for ", names(models[i+1]), " over ", names(models[i]), sep="" ))
	}
dev.off()


png("carnivora_pars_mft.png", width=800, height=400)
	plot(par_boots)
	legend("topright", levels(out$regimes), lty=c(1,2), lwd=2) 
dev.off()


