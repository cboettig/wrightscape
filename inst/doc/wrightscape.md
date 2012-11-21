% Detecting a release of constraint in Labrid fish
% Carl Boettiger, Jeremy Beaulieu and Peter C. Wainwright







# Introduction

## The release of constraint hypothesis

## Key innovations in parrotfish
[@Price2010a]


![The phylogenetic tree of Labrid fish used in this study.  We divide the tree into three clades: Wrasses, Parrotfish with an intramandibular joint and pharyngeal joint, and parrotfish that lack the intramandibular joint. Diagrams of the jaw structure for representative species from each group are shown adjacent.](figure/labrid_phylo.pdf)


## Models for a release of constraint 

The Ornstein-Uhlenbeck (OU) process is a stochastic, mean-reverting process
commonly used in the comparative phylogenetics context to model the evolution
of a trait under constraint [@Hansen1996].  The model has been generalized to 
the case of multiple optima in [@Butler2004] and recently to the case of 
differing strengths of selection [@Beaulieu2012].  This latest extension 
is perhaps most interesting when used to detect a change in selection strength
in a sub-clade following a potential innovation. This would predict and increase
in the rate disparity increases in some focal traits in the sub-clade relative
to the base rate observed.  We will refer to this model in which the strength
of stabilizing selection decreases as the release of constraint model.

Unfortunately, this basic pattern corresponds to a similar scenario that does
not involve a release of constraint.  It is commonly postulated that a trait's
evolutionary pattern may correspond best to a Brownian motion (BM) process
without no central tendency imposed by the constraint in the OU model.  If the
basic Brownian rate parameter increases at the time of the innovation, a pattern
of increased growth in disparity can still be observed without the corresponding
mechanism of a release of constraint.  Modeling a change in a Brownian rate 
parameter was first introduced in @OMeara2006 in the software *Brownie*; hence 
we will refer to this as the Brownie model. 

Though the processes are not identical, the exhibit remarkably similar patterns.
Figure 2 illustrates the Brownie and release of constraint models through 500 
replicate simulations of each.  In the top panel, BM and OU processes with 
comparable parameters are shown for reference.  In both the Brownie and release
of constraint models, the same variance has been reached at the time of the 
shift and at the time the simulation ends. The distinguishing feature in the
release of constraint model is similar to the distinction between a BM and OU
processes -- before the shift occurs, the trait dynamics have begun to approach
an equilibrium that balances the diversification process against the constraint.
This corresponds to traits values more closely reflecting a match to their 
environment than to their evolutionary history.  After the shift in selection 
occurs, the traits begin to explore outside the range previously possible under
the strong constraint.






![](figure/figure2.pdf)


By contrast, the Brownie model shows a sharper transition at this boundary. A 
close look at the figure shows a shock front at the moment of this transition.  
Whereas in the release of constraint model, the trajectories only loose their
central bias, but otherwise continue to make similar step sizes, in the Brownie
model the entire tempo of the evolutionary process has changed, taking bigger
steps in both directions. It is these subtle differences in the patterns of the
evolutionary processes implied by the different models that we seek to tease apart.
To obtain the most powerful statistical comparison between the two models that
accounts for the uncertainty in the model estimate, we use the method described
in @Boettiger2012 which uses a bootstrap simulation approach to compare likelihood
ratios of the models.  








# Results





## Estimating the release of constraint model

![](figure/figure3.pdf)


## Comparing models

![](figure/figure4.pdf)


(Show that Brownie model rejects the OUCH model/thetas?)

# Discussion

Our analysis identified two functional traits in the labrid jaw morphology that show clear evidence of a release of constraint.  


Functional innovations...



# Acknowledgements

This work was supported by a Computational Sciences Graduate Fellowship from the
Department of Energy under grant number DE-FG02-97ER25308 to CB and NSF grant
DEB-1061981 to PCW.  

# References


