# unit_test3.R

rm(list=ls()) 
require(devtools)
load_all("..")
require(ouch)
data(bimac)
tree <- with(bimac, ouchtree(node, ancestor, time/max(time), species))
data   <- as.double(
	   c(0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
            0.000000, 0.000000, 2.602690, 2.660260, 2.660260,
            2.653242, 2.674149, 2.701361, 7.161247, 13.299534,
            3.328627, 3.353407, 3.360375, 3.049273, 2.906901,
            2.980619, 2.933857, 3.975530, -3.104587, 3.346389,
            12.928524, -12.939162, 7.990720, 8.058707, -5.068053))

ws <- wrightscape(data, tree, bimac[["OU.LP"]], alpha=2.6, theta=3., sigma=.2, Xo=3.0)
lca <- lca_calc(tree)
equal <- multiOU_lik_lca(data, tree, bimac[["OU.LP"]], alpha=c(2.6,2.6,2.6), sigma=.2, theta=3, Xo=3, lca)
bad <- multiOU_lik_lca(data, tree, bimac[["OU.LP"]], alpha=c(2.6,5,2.6), sigma=.2, theta=3, Xo=3, lca)
good <- multiOU_lik_lca(data, tree, bimac[["OU.LP"]], alpha=c(2.6,1e-7,2.6), sigma=.2, theta=3, Xo=3, lca)

good > bad



	Xo <- as.double(3.0)
	alpha <- as.double(c(2.6, 2.6, 2.6))
	theta <- as.double(c(3., 3.,3))
	sigma <- as.double(c(.2,  .2, .2))
	regimes <- bimac[["OU.LP"]]

	data[is.na(data)] = 0 

	ancestor <- as.numeric(tree@ancestors)
	ancestor[is.na(ancestor)] = 0 
	ancestor <- ancestor-1  # C-style indexing

	## ouch gives cumulative time, not branch-length!!
	anc <- as.integer(tree@ancestors[!is.na(tree@ancestors)])
	lengths <- c(0, tree@times[!is.na(tree@ancestors)] - tree@times[anc] )
	branch_length <- lengths/max(tree@times)
	
	n_nodes <- length(branch_length)
	n_regimes <- length(levels(regimes))

	levels(regimes) <- 1:n_regimes
	regimes <- as.integer(regimes)-1  # convert to C-style indexing

	use_siman=0

	o<- .C("fit_model",
		as.double(Xo),
		as.double(alpha),
		as.double(theta),
		as.double(sigma),
		as.integer(regimes),
		as.integer(ancestor),
		as.double(branch_length),
		as.double(data),
		as.integer(n_nodes),
		as.integer(n_regimes),
		double(1),
		as.integer(0) 
	  )




#	#	o<- .C("calc_lik", Xo, alpha, theta, sigma, as.integer(regimes),
#	        as.integer(ancestor), branch_length, data, as.integer(n_nodes),
#		as.integer(lca), double(1)) 

