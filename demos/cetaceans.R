require(wrightscape)
require(geiger)
source("/home/cboettig/Documents/ucdavis/research/phylotrees/code/Comparative-Phylogenetics/R/data_formats.R")
path = "/home/cboettig/Documents/ucdavis/research/phylotrees/data/artiodactyla/"

tree <- read.nexus(paste(path, "artiodactyla_tree.nex", sep=""))
data <-read.csv(paste(path,"artiodactyla_mass.txt", sep=""),sep="\t")
rownames(data) <- data[,2]
matched <- treedata(tree, data)


### ape2ouch_all should really be able to handle data.frames!  Need to modify Hipp's code...
### So instead of the below, we need to do these seperately and hope it works!
#ouch_format <- ape2ouch_all(matched$phy, matched$data)

# 
mass <- log(matched$data[,3])
names(mass) <- rownames(matched$data)
ouch_format <- ape2ouch_all(matched$phy, mass)
#names(ouch_format$data) <- rownames(ouch_format$data)


regimes <-matched$data[,1]
names(regimes) <- rownames(matched$data)
ouch_regimes <- as.character( ape2ouch_all(matched$phy, regimes)$data )
ouch_regimes[is.na(ouch_regimes)] <- "anc"
ouch_regimes <- as.factor(ouch_regimes)

cetaceans <- as.integer( mrcaOUCH(ouch_format$tree@nodelabels[ ouch_regimes == levels(ouch_regimes)[3] ], ouch_format$tree) )
regimes <- paintBranches(cetaceans, ouch_format$tree)
levels(regimes) = levels(ouch_regimes)[2:3]
plot(ouch_format$tree, regime=regimes)


single_regime <- as.factor(character(ouch_format$tree@nnodes))
names(single_regime) <- ouch_format$tree@nodes

bm <- brown(ouch_format$data, ouch_format$tree)
ou1 <- hansen(ouch_format$data, ouch_format$tree, regime=single_regime, 1, 1)
ws1 <- wrightscape(ouch_format$data, ouch_format$tree, regime=single_regime, 1, 1)
ou2 <- hansen(ouch_format$data, ouch_format$tree, regime=regimes, 1, 1)
ws2 <- wrightscape(ouch_format$data, ouch_format$tree, regime=regimes, (ou2@sqrt.alpha)^2, ou2@sigma)



par_boots <- fast_boot(ws2, nboot = 200)


png("gape_pars.png", width=800, height=400)
plot(par_boots)
legend("topright", c("Wrasses", "Parrotfish"), lty=c(1,2), lwd=2) 
dev.off()






labrid_models <- list(bm = bm, ws1=ws1, ou2=ou2, ws2=ws2)

	LR <- choose_model(labrid_models, 100)
	png("labrid_protrusion.png",width=2000, height=600) 
	par(mfrow=c(1,3))
	pretty_plot(LR[[1]], main="support for OU over BM")
	pretty_plot(LR[[2]], main="support for 2 peaks over 1")
	pretty_plot(LR[[3]], main="support for differential selective strength")
	dev.off()


