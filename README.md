


# wrightscape

A fast C code implementation of maximum likelihood calculations for a multi-peak OU process with independent alpha, sigma, and theta parameters for each peak on phylogenetic trees (an extension of the ouch model).  Includes R package wrapper as interface to C code. Likelihood calculation presented in http://www.carlboettiger.info/files/phylo-covariance.pdf, later the basis for: _Jeremy M. Beaulieu, Dwueng-Chwuan Jhwueng, Carl Boettiger and Brian O’Meara, (2012). Modeling Stabilizing Selection: Expanding the Ornstein-Uhlenbeck Model of Adaptive Evolution, Evolution 66 (8) 2369-2383. doi:10.1111/j.1558-5646.2012.01619.x_


## File overview

`regimes.c`: 
	Test the regimes data structure
	Modify to handle a single regime
	Implement a transition density calculation where transitions occur at nodes alone  <<< 
	Pagel & Meade, Huelsenbeck papers -- what about transitions at nodes?  
	
`linear.c`: 
	BM not really subset of OU, as alpha = 0 not acceptible OU calculation.
	Modify the OU model to handle appropriately small alpha with Taylor Expansion or as BM

`max likelihood`:
	Implement and test maximum likelihood searches over parameters against existing methods

`mcmc.c`:
	Start implementing a basic MCMC solver!


`matrix_method.c`:
	nonlinear transition densities?  
	Test the stationary distribution?



### Error handling  

This should probably be done mostly at the R level, enforcing the strict formats and conventions of the input (i.e. if not doing regimes, all states should be listed at 0.  States must be numbered incrementally from 0 to n-1, where n is the number of regimes present.)





#### STEPPING STONES ####

- Recursive Joint Prob method using BM, compare to analyticals
- Bayesian OU vs BM in RJMCMC
- Use Pagel's BayesTraits run to "paint tree" with uncertainty, then repeat OUCH analysis with uncertainty in branch color



