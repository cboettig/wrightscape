# labrid example

# This data has not been released
require(wrightscape)
require(geiger)
source("/home/cboettig/Documents/ucdavis/research/phylotrees/code/Comparative-Phylogenetics/R/data_formats.R")

#path = "/home/cboettig/Documents/research/phylotrees/data/labrids/"
path = "../../data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
labrid_data <-read.csv(paste(path,"labrid.csv", sep=""))
rownames(labrid_data) <- labrid_data[,1]

# We'll just use the fin angle data
fin_angle <- labrid_data$FinAngle
names(fin_angle) <- labrid_data[,1]

# Convert the data, dropping unmatched tips
labrid <- treedata(labrid_tree, fin_angle)
names(labrid$data) <- rownames(labrid$data)
labrid <- ape2ouch_all(labrid$phy, labrid$data)

# We'll paint the parrotfish a different regime
parrotfish <- labrid$tree@nodelabels[c(pmatch("Scarus_spinus", labrid$tree@nodelabels), pmatch("Crypto", labrid$tree@nodelabels) )]
parrotfish_clade <- as.integer(mrcaOUCH(parrotfish, labrid$tree))

parrotfish_regime <- paintBranches(parrotfish_clade, labrid$tree)

single_regime <- as.factor(character(labrid$tree@nnodes))
names(single_regime) <- labrid$tree@nodes

bm <- brown(labrid$data, labrid$tree)
ou1 <- hansen(labrid$data, labrid$tree, regime=single_regime, 1, 1)
ws1 <- wrightscape(labrid$data, labrid$tree, regime=single_regime, 1, 1)
ou2 <- hansen(labrid$data, labrid$tree, regime=parrotfish_regime, 1, 1)
ws2 <- wrightscape(labrid$data, labrid$tree, regime=parrotfish_regime, (ou2@sqrt.alpha)^2, ou2@sigma)

labrid_models <- list(bm = bm, ws1=ws1, ou2=ou2, ws2=ws2)

	LR <- choose_model(labrid_models, 100)
	png("wrightscape_labrid.png",width=2000, height=600) 
	par(mfrow=c(1,3))
	pretty_plot(LR[[1]], main="support for OU over BM")
	pretty_plot(LR[[2]], main="support for 2 peaks over 1")
	pretty_plot(LR[[3]], main="support for differential selective strength")
	dev.off()


