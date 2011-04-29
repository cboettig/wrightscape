# labrid example

require(wrightscape)
require(pmc)
require(socialR)
tag="phylogenetics wrightscape labrids"

# This data has not been released
path = "../data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
#fin_data <-read.csv(paste(path,"labrid.csv", sep=""))
diet_data <- read.csv(paste(path,"labriddata_parrotfish.csv", sep=""))

corrected_data <- diet_data

# Use the simple names from Price et al 2010
traitnames <- c("Species", "group", "gape", "prot", "bodymass", "AM", "SH", "LP", "close", "open", "kt")
names(corrected_data) <- traitnames

# Lengths are log transformed 
corrected_data[["gape"]] <- log(corrected_data[["gape"]])
corrected_data[["prot"]] <- log(corrected_data[["prot"]])
#masses are log(cube-root) transformed
corrected_data[["bodymass"]] <- log(corrected_data[["bodymass"]])/3
corrected_data[["AM"]] <- log(corrected_data[["AM"]])/3
corrected_data[["SH"]] <- log(corrected_data[["SH"]])/3
corrected_data[["LP"]] <- log(corrected_data[["LP"]])/3
# ratios are fine as they are


## We'll perform Revell's regression on all size/mass (non-ratio) traits
# 
size_col <- 5
cols_to_size_correct <- c(3:4,6:8)
ape <- treedata(labrid_tree, corrected_data, corrected_data[,1])
# Annoying formating to get a labeled matrix with numerics, not characters
Y <- matrix(as.numeric(ape$data[,cols_to_size_correct]), ncol=length(cols_to_size_correct))
colnames(Y) <- colnames(ape$data[,cols_to_size_correct])
rownames(Y) <- ape$data[,1] 
x <-  Y[,size_col]
names(x) <- rownames(Y)

require(RevellExtensions)
out <- phyl.resid(labrid_ape$phy,x, Y)
## phyl.resid changes order of species listing!

#Could duplicate one column and change the others, and merge by that.  OR just use by="row.names"
resid <- data.frame(out$resid, rownames(out$resid))
colnames(resid) <- c("gape_cor", "prot_cor", "AM_cor", "SH_cor", "LP_cor", "Species")
test <- merge(ape$data, resid)
# The above is a good double-check on the order 

traits <- merge(ape$data, out$resid, by="row.names")
# columns that are transformed now have gape.x for untransformed, gape.y for transformed.  

## Merge sanity test
#head(test)
#head(traits)



# format data gets regimes from column specified in "regimes" (e.g. this is a column id, not a # of regimes)
# This also converts the tree and data into ouch format

## Could just hand it all traits, but these are just the tranformed and size-corrected ones
labrid <- format_data(labrid_tree, traits[c(6,10:17)], species_names=traits[,1],  regimes = 3)  
names(labrid$data)


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



cpu <- 9
nboot <- 160
sfInit(parallel=TRUE, cpu=cpu)
sfExportAll()
sfLibrary(wrightscape)  # need all this just to export wrightscape?
sfLibrary(pmc)

# try only for these, or try for all traits
#c(3,4,5,9,10,11)

sfSapply(1:9, function(i){
	trait_name <- names(labrid$data)[i]	
	trait <- labrid$data[i]

	bm <- brown(trait, labrid$tree)
	ou1 <- hansen(trait, labrid$tree, regime=labrid$noregimes, .01, .01)
	ou2_phar <- hansen(trait, labrid$tree, regime=labrid$regimes, .01, .01)
  ou2_intra <- hansen(trait, labrid$tree, regime=intramandibular, .01, .01)
  ou3 <- hansen(trait, labrid$tree, regime=two_shifts, ou2_phar@sqrt.alpha, sigma=ou2_phar@sigma )

  ## Note that this will converge poorly with the .01, .01 starting conditions
  #  ou3 <- hansen(trait, labrid$tree, regime=two_shifts, .01, .01 )
  #  loglik(ou3) - loglik(ou2_phar)



  ## Compare intramandibular and pharyngeal joint paintings
  ##  Note that we use the more restricted estimates from ouch to seed the search

#  pharyngeal vs intramadibular regimes

  ouch_phar <- ouch(trait, labrid$tree, regime=labrid$regimes, alpha=(ou2_phar@sqrt.alpha)^2, sigma=ou2_phar@sigma)
	brownie_phar <- brownie(trait, labrid$tree, regime=labrid$regimes, sigma=ou2_phar@sigma)
	wright_phar <- wright(trait, labrid$tree, regime=labrid$regimes, alpha=(ou2_phar@sqrt.alpha)^2, sigma=ou2_phar@sigma)
  release_phar <- release_constraint(trait, labrid$tree, regime=labrid$regimes, alpha=(ou2_phar@sqrt.alpha)^2, sigma=ou2_phar@sigma)

  ouch_intra <- ouch(trait, labrid$tree, regime=intramandibular, alpha=(ou2_phar@sqrt.alpha)^2, sigma=ou2_phar@sigma)
	brownie_intra <- brownie(trait, labrid$tree, regime=intramandibular, sigma=ou2_intra@sigma)
	wright_intra <- wright(trait, labrid$tree, regime=intramandibular, alpha=(ou2_intra@sqrt.alpha)^2, sigma=ou2_intra@sigma)
  release_intra <- release_constraint(trait, labrid$tree, regime=intramandibular, alpha=(ou2_intra@sqrt.alpha)^2, sigma=ou2_intra@sigma)

  ouch_twoshifts <- ouch(trait, labrid$tree, regime=two_shifts, alpha=(ou3@sqrt.alpha)^2, sigma=ou3@sigma)
	brownie_twoshifts <- brownie(trait, labrid$tree, regime=two_shifts, sigma=ou3@sigma)
	wright_twoshifts <- wright(trait, labrid$tree, regime=two_shifts, alpha=(ou3@sqrt.alpha)^2, sigma=ou3@sigma)
  release_twoshifts <- release_constraint(trait, labrid$tree, regime=two_shifts, alpha=(ou3@sqrt.alpha)^2, sigma=ou3@sigma)
 
  loglik(wright_twoshifts)-loglik(wright_intra)


  results <- matrix(NA, nrow=4, ncol=3, dimnames = list(c("ouch", "brownie", "release", "wright"), c("phar", "intra", "twoshifts")))

  results[1,1] <- loglik(ouch_phar)
	results[2,1] <- loglik(brownie_phar)
  results[3,1] <- loglik(release_phar)
	results[4,1] <- loglik(wright_phar)

  results[1,2] <- loglik(ouch_intra )
	results[2,2] <- loglik(brownie_intra )
  results[3,2] <- loglik(release_intra)
	results[4,2] <- loglik(wright_intra )

  results[1,3] <- loglik(ouch_twoshifts )
	results[2,3] <- loglik(brownie_twoshifts)
  results[3,3] <- loglik(release_twoshifts)
	results[4,3] <- loglik(wright_twoshifts)

#  barplot(results, xlim=c(0, 80), col=c("thistle", "khaki", "pink","palegreen"), horiz=TRUE, beside=TRUE)
  social_plot(barplot(t(results), xlim=c(0, 80), col=c("thistle", "khaki", "palegreen"), horiz=TRUE, beside=TRUE, main=trait_name), tag=tag, comment=trait_name)

}




 # Can we do better?
# wright_twoshifts_ <- wright(trait, labrid$tree, regime=two_shifts, alpha=c(wright_intra$alpha, wright_intra$alpha[2]), sigma=c(wright_intra$sigma, wright_intra$sigma[2]))


## Doing the likelihood optimization in C instead.  needs robustness testing, convergence conditions still
#  ws2_phar <- wrightscape(trait, labrid$tree, regime=labrid$regimes, (ou2_phar@sqrt.alpha)^2, ou2_phar@sigma, theta=ou2_phar@theta[[1]])
#  ws2_intra <- wrightscape(trait, labrid$tree, regime=intramandibular, (ou2_intra@sqrt.alpha)^2, ou2_intra@sigma, theta=ou2_intra@theta[[1]])
#  ws2_twoshifts <- wrightscape(trait, labrid$tree, regime=intramandibular, (ou3@sqrt.alpha)^2, ou3@sigma, theta=ou3@theta[[1]])


## Make this pretty using pmc and ape plot tools
plt <- function(){
  phylo_phar <- convert(ou2_phar)
  phylo_intra <- convert(ou2_intra)
  phylo_twoshifts <- convert(ou3)
  ou3 <- hansen(trait, labrid$tree, regime=two_shifts, ou2_phar@sqrt.alpha, sigma=ou2_phar@sigma )
  par(mfrow=c(1,3))
  plot(phylo_phar, edge.color = treepalette(phylo_phar), edge.width=2, cex=1.0, show.tip.label=FALSE)
  plot(phylo_intra, edge.color = treepalette(phylo_intra), edge.width=2, cex=1.0, show.tip.label=FALSE)
  plot(phylo_twoshifts, edge.color = treepalette(phylo_twoshifts), edge.width=2, cex=1.0, show.tip.label=FALSE)
}
#social_plot(plt(), file="labrids.png", tag=tag)


#plot(labrid$tree, regimes=two_shifts)


## Do some modelchoice
#out <- montecarlotest(brownie_phar, release_phar, cpu=cpu,nboot=nboot, GetPar=F) 
#social_plot(plot(out), tag="phylogenetics wrightscape labrids", comment="brownie vs release on pharyngeal shift pt, trait = gape")
#out2 <- montecarlotest(release_phar, release_intra, cpu=cpu,nboot=nboot, GetPar=F) 
#social_plot(plot(out2), tag="phylogenetics wrightscape labrids", comment="release on pharyngeal vs intramandibular shift, trait=gape")
#out3 <- montecarlotest(release_intra, release_twoshifts, cpu=cpu,nboot=nboot, GetPar=F) 
#social_plot(plot(out3), tag="phylogenetics wrightscape labrids", comment="release on pharyngeal vs intramandibular shift, trait=gape")



#})



