#
# Adapted from "maticce" by Andrew Hipp, now obsolete

# ------------------------------------------
# FUNCTIONS FOR TRAVERSING AN S4 OUCH TREE #
# ------------------------------------------

isMonophyletic <- function(tree, taxa) {
# returns T or F on whether a group of taxa is monophyletic in an ouch tree
  if(length(taxa) == 1) return(taxa %in% tree@nodelabels[tree@term])
  else(return(identical(sort(taxa), sort(nodeDescendents(tree, mrcaOUCH(taxa, tree))))))
}

nodeDescendents <- function(tree, startNode) {
## Recursive function to find all the descendents of a node on an 'ouchtree' object
  startNode <- as.character(startNode) # just to be safe
  daughterBranches <- as.character(tree@nodes[tree@ancestors %in% startNode])
  nodeNames <- tree@nodelabels[tree@nodes %in% daughterBranches]
  if(!identical(as.character(daughterBranches), character(0))) {
    for(i in daughterBranches) nodeNames <- c(nodeNames, nodeDescendents(tree, i))
  }
  return(nodeNames[nodeNames %in% tree@nodelabels[tree@term]])
}
  
mrcaOUCH <-
# Finds most recent common ancestor for a vector of tips by:
#  1. Creating a vector of ancestral nodes leading to each tip
#  2. Creating an intersection set of ancestral nodes for all taxa by intersecting taxa successively with the last intersection set
#  3. Returning the node of the final intersection set that has the highest time
# Arguments:
#  "node" "ancestor" "times" "species" = the standard tree specification vectors of the OUCH-style tree
#  "cladeVector" = vector of species for which you want to find the most recent common ancestor
# Value: the node number of the most recent common ancestor
function(cladeVector, tree) {
  ## ------------------ begin ouchtree block -----------------
  ## check to see if tree inherits 'ouchtree'
  if (!is(tree,'ouchtree')) 
	stop(paste('This function has been rewritten to use the new S4 ', sQuote('ouchtree'), ' class.',
	'\nYou can generate a tree of this class by calling ', sQuote('ouchtree()'), '.', sep = ""))
  ## get the vectors we need:
  ancestor <- tree@ancestors # class = "character"
  node <- tree@nodes # class = "character"
  species <- tree@nodelabels # class = "character" -- note that nodelabels is more general than this indicates and the name should be changed throughout at some point
  times <- tree@times # class = "numeric"
  ## ------------------ end ouchtree block -------------------
  
  if(length(cladeVector) == 1) return(tree@nodes[tree@nodelabels == cladeVector])
  else {
    tips = match(cladeVector, species) 
    listOfAncestorLines = lapply(tips, ancestorLine, tree = tree) # 10 nov 08: this is identical to the appropriate subset of tree@lineages
    latestMatch = listOfAncestorLines[[1]]
    for (i in listOfAncestorLines) {
      latestMatch = i[match(latestMatch, i, nomatch = 0)] }
    timesVector = times[as.integer(latestMatch)]
    if(length(timesVector) == 1) {
      if (is.na(timesVector)) mrca = "1"
        else mrca = timesVector}
      else mrca = latestMatch[match(max(as.double(timesVector), na.rm = TRUE), timesVector)]
    return(mrca) 
  }
}

ancestorLine <-
# Creates a vector of ancestral nodes for a tip
# Arguments:
#  "tree" = an ouch-style (S4) tree
#  "tip" = the tip node to trace back
# Value: a vector of nodes leading from a tip to the root
# 10 nov 08: changed to just grab tree@lineages and make a vector that fits the old code
function(tip, tree) {
  ## check to see if tree inherits 'ouchtree'
  if (!is(tree,'ouchtree')) 
	stop(paste('This function has been rewritten to use the new S4 ', sQuote('ouchtree'), ' class.',
	'\nYou can generate a tree of this class by calling ', sQuote('ouchtree()'), '.', sep = ""))
  tip <- as.numeric(tip)
  nodesVector <- c(as.character(tree@lineages[[tip]][2:length(tree@lineages[[tip]])]), NA)
  return(nodesVector) 
  }
