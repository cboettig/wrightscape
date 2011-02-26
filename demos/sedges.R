# Example using Andrew Hipp's sedges dataset, from the maticce package.  

	source("/home/cboettig/Documents/ucdavis/research/phylotrees/code/Comparative-Phylogenetics/R/data_formats.R")
	require(wrightscape)
	require(maticce)
	data(carex)
	attach(carex)

	treedata <- ape2ouch_all(ovales.tree, ovales.data)
	ou2_regimes <- paintBranches(list(ovales.nodes[[2]]), treedata$tree)
	ou1_regimes <- as.factor(rep(1, length(ou2_regimes) ))
	names(ou1_regimes) <- names(ou2_regimes)

# Fit the models
	bm <- brown(treedata$data, treedata$tree)
	ws1 <- wrightscape(treedata$data, treedata$tree, ou1_regimes, 1, 1)
	ou2 <- hansen(treedata$data, treedata$tree, ou2_regimes, 1, 1)
	ws2 <- wrightscape(treedata$data, treedata$tree, ou2_regimes, (ou2@sqrt.alpha)^2, ou2@sigma)

	model_list <- list(bm = bm, ws1 = ws1, ou2 = ou2, ws2 = ws2)
	LR <- choose_model(model_list, 100)
	par(mfrow=c(1,3))
	pretty_plot(LR[[1]], main="support for OU over BM")
	pretty_plot(LR[[2]], main="support for 2 peaks over 1")
	pretty_plot(LR[[3]], main="support for differential selective strength")


