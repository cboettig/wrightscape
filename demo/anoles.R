# Example using the Anoles of the Lesser Antilles from Butler & King 2004

	require(wrightscape)
	
# Load the data
	data(bimac)
	tree <- with(bimac,ouchtree(nodes=node,ancestors=ancestor,times=time/max(time),labels=species))

# Fit the models
	bm <- brown(log(bimac['size']), tree)
	ou1 <- hansen(log(bimac['size']), tree, bimac['OU.1'], 1, 1)
	ws1 <- wrightscape(log(bimac['size']), tree, bimac['OU.1'], 1, 1)
	ou2 <- hansen(log(bimac['size']), tree, bimac['OU.LP'], 1, 1)
	ws2 <- wrightscape(log(bimac['size']), tree, bimac['OU.LP'], (ou2@sqrt.alpha)^2, ou2@sigma)

# Bootstrap comparisons -- slow! 
	model_list <- list(bm = bm, ws1 = ws1, ou2 = ou2, ws2 = ws2)
	LR <- choose_model(model_list, 100)

# Plot results of likelihood ratio test
	par(mfrow=c(1,3))
	pretty_plot(LR[[1]], main="support for OU over BM")
	pretty_plot(LR[[2]], main="support for multiple peaks over 1")
	pretty_plot(LR[[3]], main="support for differential selective strength")
	
# some "by hand" methods available
	LR <- LR_bootstrap(ou2, ws2, n=100)
	plot(LR)
	x <- simulate(ws2)
	up_w <- update(ws2, x$rep.1)
	up_h <- update(ou2, x$rep.1)



