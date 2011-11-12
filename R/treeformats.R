# treeformats.R
# Convert toggles between ouch and ape format trees 

convert <- function(ot, regimes=NULL){
	if(is(ot, "ouchtree")){

		n <- ot@nnodes
		# find the tips vs internals
		anc <- as.integer(ot@ancestors[!is.na(ot@ancestors)])
		internal <- sort(anc)[seq(2,n-1, by=2)]
		tmp <- integer(n)
		tmp[internal] = 1
		tips <- which(tmp == 0)

		root <- which(is.na(ot@ancestors))
		internal <- internal[internal!=root] #remove root form internal list
		
		new_ids <- integer(n)
		new_ids[tips] <- 1:length(tips)
		new_ids[root] <- length(tips)+1
		new_ids[internal] <- (length(tips)+2):n
		
		new_ancestor <- new_ids[as.integer(ot@ancestors)]

		edge <- matrix(NA, n-1, 2)
		edge[,1] <- new_ancestor[!is.na(new_ancestor)]
		edge[,2] <- new_ids[!is.na(new_ancestor)]

		anc <- as.integer(ot@ancestors[!is.na(ot@ancestors)])
		lengths <- ot@times[!is.na(ot@ancestors)] - ot@times[anc]
	
		labels <- ot@nodelabels[tips]

		tree <- list(edge=edge, Nnode = (n-1)/2, tip.label = labels, edge.length= lengths )
		class(tree) <- "phylo"

		if (is(ot, "hansentree")) {
			regimes <- ot@regimes[[1]][-1]
			tree$regimes <- regimes
		}

    if(!is.null(regimes)){
      tree$regimes <- regimes[-1]
    }

	} else	if (is(ot, "phylo")){ 
			tree <- ape2ouch(ot)
		}
#	plot(tree)
	tree
}




# coloring for trees 
treepalette <- function(apetree, colormap = c("rainbow", "heat.colors", "terrain.colors", "topo.colors", "cm.colors", "gray"), custom=NULL, rev=FALSE ){ 
	colormap <- match.arg(colormap)
	if(colormap=="rainbow") 
    levels(apetree$regimes) <- rainbow(length(levels(apetree$regimes)))
	if(colormap=="heat.colors") 
    levels(apetree$regimes) <- heat.colors(length(levels(apetree$regimes)))
	if(colormap=="terrain.colors") 
    levels(apetree$regimes) <- terrain.colors(length(levels(apetree$regimes)))
	if(colormap=="topo.colors") 
    levels(apetree$regimes) <- topo.colors(length(levels(apetree$regimes)))
	if(colormap=="cm.colors") 
    levels(apetree$regimes) <- cm.colors(length(levels(apetree$regimes)))
	if(colormap=="gray"){
    n <-length(levels(apetree$regimes))+2
    levels(apetree$regimes) <- gray((1:n)/n)
  }
  if(!is.null(custom))
    levels(apetree$regimes) <- custom
	if(rev==TRUE) 
    levels(apetree$regimes) <- rev(levels(apetree$regimes))
	as.character(apetree$regimes)
} 



