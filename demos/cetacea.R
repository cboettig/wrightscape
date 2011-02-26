# cetacea.R   
# Body size and diet.  No artiodactyle set included 
# Not sure that the regimes specified determine monophyletic groups


nboot <- 100
cpu = 8

dataset <- "cetacea"

require(wrightscape)
#path = paste("~/data/", dataset, sep="")
path  = paste("../../../data/", dataset, sep="")
path = paste("/home/cboettig/Documents/ucdavis/research/phylotrees/data/", dataset, sep="")



tree <- read.nexus(paste(path, "/", dataset, ".nex", sep=""))
data <- read.csv(paste(path, "/", dataset, ".txt", sep=""), sep="\t")
out <- format_data(tree, data, species_names = data[,1], regimes = data[,3])

size <- log(out$data[,2])
names(size) <- rownames(out$data)

bm <- brown(size, out$tree )
ou1 <- hansen(size, out$tree, regimes=out$noregimes, 1,1 )
ou2 <- hansen(size, out$tree, regimes=out$regimes, 1,1 )
ws2 <- wrightscape(size, out$tree, regime=out$regimes, (ou2@sqrt.alpha)^2, ou2@sigma, theta=ou2@theta[[1]])

models <- list(bm = bm, ou1 = ou1, ou2 = ou2, ws2 = ws2)
LR <- choose_model(models, nboot=nboot, cpu=cpu)
	pdf(paste(dataset, "_models.pdf", sep=""),width=1500, height=600) 
	par(mfrow=c(1,(length(models)-1) ))
	for(i in 1:(length(models)-1) ){
		pretty_plot(LR[[i]], main=paste("support for ", names(models[i+1]), " over ", names(models[i]), sep="" ))
	}
dev.off()


par_boots <- fast_boot(ws2, nboot = nboot, cpu=cpu)
	pdf(paste(dataset, "_pars.pdf",sep=""), width=800, height=400)
	plot(par_boots)
	legend("topright", levels(out$regimes), lty=c(1,2), lwd=2) 
dev.off()


