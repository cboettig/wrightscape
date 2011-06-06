# paint_labrid_regimes.R

## PAINTING REGIMES: Having taken care of traits, we paint on the 3 regime models.  
# Select common ancestor of a Chlorurus and a Hipposcarus as the changepoint
intra_ancestor <- mrcaOUCH(c("Chlorurus_sordidus", "Hipposcarus_longiceps"), labrid$tree)
intramandibular <- paintBranches(intra_ancestor, labrid$tree, c("other","intramandibular"))

# Select common ancestor for all parrot fish
pharyngeal_ancestor<-mrcaOUCH(c("Bolbometopon_muricatum", "Sparisoma_radians"), labrid$tree)
pharyngeal <- paintBranches(pharyngeal_ancestor, labrid$tree, c("other","pharyngeal"))

## A third painting involves shifts at both points
two_shifts <- paintBranches(c(pharyngeal_ancestor, intra_ancestor), labrid$tree, c("wrasses", "pharyngeal", "intramandibular") )

## This leaves the branch on which the second transition occurs unspecified (fourth regime).  
## We have to fix this manually
two_shifts[as.numeric(intra_ancestor)] <- "intramandibular"
two_shifts <- as.factor(as.character(two_shifts))
names(two_shifts) <- names(intramandibular)


