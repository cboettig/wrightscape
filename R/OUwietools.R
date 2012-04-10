
#' returns a phylogenetic tree painted for ouwie and data
#' 
#' @param phy a phylogenetic tree of class phylo (ape)
#' @param clades: a list of clades defined by species pairs,
#' e.g. list(c(sp1, sp2), c(sp3, sp4)), or just a vector c(sp1, sp2)
#' the descendents of the MRCA of sp1 and sp2 define the regime
#' later clades are painted over earlier ones, so should be in 
#' order of nesting or otherwise non-overlapping. 
#' @param data the data.frame of characters (see data in ?treedata)
#' @param show_plot logical, should I plot the resulting tree
#' @return a list with the phy, updated to have regimes at node labels
#' (for OUwie) and the data matrix formatted with regimes, for OUwie,
#' and a list of colors corresponding to edges, that can be used for plotting.  
#' @examples
#' data(labrids)
#' myphy <- paint_phy(ape.phy, ape.dat, 
#'                    list(c("Bolbometopon_muricatum", "Sparisoma_radians"), 
#'                    c("Chlorurus_sordidus", "Hipposcarus_longiceps")))
#' plot(myphy$phy, edge.color=myphy$colors, type="fan", 
#' show.tip.label=FALSE, edge.width=2)
#' 
#' @import ape
#' @import geiger
#' @export 
paint_phy <- function(phy, data, clades, show_plot=TRUE){

 # drop unmatched tips & characters
 # may produce difficulties on already-dropped trees?
  pruned <- treedata(phy, data)
  phy <- pruned$phy
  data <- pruned$data

  # one regime or multiple?
  if(is.list(clades)){
    sp1 <- sapply(clades, function(x) x[1])
    sp2 <- sapply(clades, function(x) x[2])
    n <- length(clades)
  } else if(is.character(clades)){
    sp1 <- clades[1]
    sp2 <- clades[2]
    n <- 1
  } else {
    error("clades input not recognized")
  }

  desc <- vector("list", length=n)
  regimes <- rep(0, length(phy$edge[,1]))
  colors <- rep(1, length(phy$edge[,1]))

  for(i in 1:n){
    # Create the painting
    desc[[i]] <- get_clade(phy, sp1[i], sp2[i])
    regimes[desc[[i]]] <- i
    colors[desc[[i]]] <- i+1
  }

  # pick out tips and nodes
  tips <-  !(phy$edge[,2] %in% phy$edge[,1]) 
  regime_tips <- regimes[tips]
  regime_nodes <- c(0, regimes[!tips]) # adds root back on

  # stick them on the tree
  phy$node.label <- as.character(regime_nodes)

  Reg <- data.frame(Reg=as.integer(regime_tips),
                    row.names=phy$tip.label)
  data <- merge(Reg, as.data.frame(data), by="row.names")
  names(data)[1] <- "Genus_species"

  # optionally plot the tree to verify the painting
  if(show_plot){
    plot(phy, edge.color=colors)
    nodelabels(pch=21, bg=c("black", colors[!tips]))
    tiplabels(pch=21, bg=colors[tips])
  }

  list(phy=phy, data=data, colors=colors)
}


#' A function that returns the edge numbers of all branches 
#' descended from the mrca of sp1 and sp2.  
#' @param phy a phylogenetic tree of class phylo (ape)
#' @param sp1 the tip.label of the 1st species
#' @param sp1 the tip.label of the 2nd species 
#' the descendents of the MRCA of sp1 and sp2 define the regime
get_clade <- function(phy, sp1, sp2){
  M <- mrca(phy)
  who <- M[sp1, sp2]
  desc <- which(phy$edge[,1] %in% who)
  l1 <- length(desc)
  l2 <- length(desc)+1
  while(l1!=l2){
    l1 <- length(desc)
    who <- phy$edge[desc, 2]
    desc <- unique(c(desc, which(phy$edge[,1] %in% who)))
    l2 <- length(desc)
  }
  desc
}

