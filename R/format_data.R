#' format data in ape format into ouch format
#' @param tree a phylogenetic tree of class "phylo", ape format
#' @param traits a numeric with trait values, or a matrix or data frame of traits, rownames matching species or handed in
#' @param species_names in the order of entries in traits, if not given in rownames.  
#' @param regimes the column in traits containing regime labels
#' @return the ouch-formatted tree, traits, and regimes  
#' @seealso \code{\link{convert}} to toggle between formats, including regime paintings
#' @details Should become an internal function to handle data conversion to ape-type 
#' @export 
format_data <- function(tree, traits, species_names = NULL, regimes = NULL ){
	require(geiger)

# Function checks that tree and trait match and convert them into a format used by wrightscape
# Function also will code tree by finding the common ancestor of all species with matching entry specified in the regimes list and assigning that codename as the regime of all descendents of that ancestor.  May not handle conflicts if corresponding to overlapping clades.  Alternatively, the regimes can be specified directly in ouch format.   

	# Args:
		# tree can be a path to a nexus file, an object of class "phylo", or an object of class "ouchtree"


	if( is(tree, "character" ) ){ 
		tree <- read.nexus(tree) 
	} else if (is(tree, "ouchtree")) {
		# uses my ouch2ape tree conversion script
		#tree <- convert(tree)
	}
	if( !is(tree, "phylo") ) { stop("Problem with tree format") }


	# Figure out if species names is already attached to traits 
	if(is.null(species_names) ){
		if( !is.null(names(traits)) & is(traits, "numeric") ){
			species_names <- names(traits)
		} else if( !is.null(rownames(traits)) ){ 
			species_names <- rownames(traits) 
		} else { 
			stop("Species names not found")
		}
	}

	## attach species names to traits
	if(is(traits, "numeric") ){
		names(traits) <- species_names
	} else if (is(traits, "data.frame") | is(traits, "matrix")){
	rownames(traits) <- species_names
	} else { stop("traits format unrecognized") }

	# drop missing taxa using geiger's function
	matched <- treedata(tree, traits, species_names)
	# treedata returns a matrix which loses class type for trait data.frame
	# so we restore the the dataframe structure

	if( is(matched$data, "matrix") ){ 
		matched$data <- as.data.frame(matched$data, stringsAsFactors=FALSE)
		if(is(traits, "data.frame")){
		for(i in 1:length(traits)){
				tmp <- class(traits[[i]])
				# each datarow that isn't numeric should be a character string (for regime, etc)
				if(tmp == "factor") tmp <- "character"
				class(matched$data[[i]]) <- tmp
			}
		}
	} 
	traits <- matched$data
	species_names <- rownames(matched$data) 

	# Convert to OUCH format 
	tree <- ape2ouch(matched$phy) 
   
	# Makes sure data is reformatted to ouch format matching the tree
	if( is(traits, "data.frame") ){
		message("traits ouch-formatted as data.frame")
		if(length(traits)>1){
		    dataIn <- traits[match(tree@nodelabels, species_names),]
		} else {
			#hack to get around size-1 data-frames 
			traits[[2]] = NA
			dataIn <- traits[match(tree@nodelabels, species_names),]
			dataIn[[2]] <- NULL
		}
	} else if(is(traits, "numeric") ){ 
message("traits ouch-formatted as numeric")
	    dataIn <- traits[match(tree@nodelabels, species_names)]
	} else { stop(paste("data of class", class(traits), "not recognized")) }
	rownames(dataIn) <- tree@nodes

	if(is.numeric(regimes)){	
		R <- compute_regimes(tree, traits, species_names, regimes)
    nr <- R$noregimes
	} else {
		regimes= as.factor(rep(" ", length=tree@nnodes))
		names(regimes) <- tree@nodes 
		R <- list(regimes=regimes) 
    nr <- R$regimes
	}

	list(tree=tree, data=dataIn, regimes=R$regimes, noregimes=nr)
} 


compute_regimes <- function(tree, traits, species_names, regimes){
# regimes can be an integer specifying the column in traits dataframe 
# Internal function for convert_data
# takes ouch_formatted tree and traits,  

	## attach species names to regimes
	if(is.null(regimes) ){
		if( is(traits, "data.frame") | is(traits, "matrix") ){
			n <- dim(traits)[1]
		} else if (is(traits, "numeric")){
			n <- length(traits)
		}
		regimes <- as.factor(character(n))
		names(regimes) <- species_names
	} else if(length(regimes) == 1) {
		regimes <- traits[,regimes]
	}

	# order regimes into ouch format
	regimes <- regimes[match(tree@nodelabels, species_names)]

	# add an ancestral regime instead of NA, which confuses the mrcaOUCH function
	regimes <- as.character(regimes)
	regimes[is.na(regimes)] <- "anc"
	regimes <- as.factor(regimes)

	# add a noregimes case for ou1 fits
	noregimes <- as.factor(character(length(regimes)))

	# Paint regimes by clade using functions from maticce
	clades <- which(levels(regimes) != "anc")
	clade_ancestors <- integer(length(levels(regimes))-1)
	clade_names <- character(length(clade_ancestors) )
	k <- 1
	for(i in clades ){
		ancestor <-	as.integer(
			mrcaOUCH( tree@nodelabels[ regimes == levels(regimes)[i] ], tree)
		)
		clade_ancestors[k] <- ancestor
		clade_names[k] <- levels(regimes)[i]
		k <- k+1
	}
	## paintBranches needs the regime names in order
	A <- data.frame(clade_ancestors, clade_names)
	A <- A[order(A$clade_ancestors), ]
	regimes <- paintBranches(A$clade_ancestors, tree, regimeTitles=as.character(A$clade_names))

	# both regimes need node numbers as their names for OUCH
	names(noregimes) <- names(regimes)
	list(regimes=regimes, noregimes=noregimes)
}












## adapted from maticce
ape2ouch_all <-
 function(tree, characterStates){
	if( is(tree, "phylo")){
		tree <- ape2ouch(tree) 
	} else if(is(tree, "ouchtree") ) {
		# already an ouchtree format
	} else { stop(paste("Input tree type not recognized" )) }
  ## Check character states to make sure that they are either named and match names in the trees, or are the same length as the tips
    dataFlag <- NULL
    stopFlag <- FALSE
    tree 
    terminals <- tree@nodelabels[(tree@nnodes - tree@nterm + 1):tree@nnodes]
    if(any(FALSE %in% (terminals %in% names(characterStates)))) {
      message(paste("Not every terminal branch in tree has a corresponding name in", sQuote("characterStates")))
      if(length(characterStates) == tree@nterm) {
        message("Data assumed to be in the same order as terminals")
        dataFlag <- 'sameOrderTerminals' 
        }
      if(length(characterStates) == tree@nnodes) {
        message("Data assumed to be in the same order as nodes;\nany data not associated with a terminal branch will be ignored")
        dataFlag <- 'sameOrderNodes'
        }
      if(identical(dataFlag, NULL)) stopFlag <- TRUE
      message("-------------------\n")
      }
    else dataFlag <- 'named'
    if(stopFlag) stop("Correct discrepancies between trees and data and try again!")
    
    ## make sure data fits the tree
    dataIn <- NULL
    if(dataFlag == 'sameOrderTerminals') dataIn <- c(rep(NA, tree@nnodes - tree@nterm), characterStates)
    if(dataFlag == 'sameOrderNodes') dataIn <- characterStates
    if(dataFlag == 'named') dataIn <- characterStates[match(tree@nodelabels, names(characterStates))]
    else names(dataIn) <- tree@nodes

	list(tree=tree, data=dataIn)
} 



