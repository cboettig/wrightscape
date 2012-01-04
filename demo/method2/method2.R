# ouwie_parrotfish.R
# Author: Carl Boettiger <cboettig@gmail.com>
# Date: 19 October 2011
require(OUwie)


data(tworegime)

## cannot run this with nodes
tree$node.label[1:length(tree$node.label)] <- 1
#tree$node.label[length(tree$node.label)] <- 2

#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label)) 
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
plot(tree) 
nodelabels(pch=21, bg=select.reg)

#To see the first 5 lines of the data matrix to see what how to
#structure the data:
trait[,"X"] <- rnorm(dim(trait)[1])

require(geiger)
tree <- lambdaTree(tree, 0.0001)
plot(tree)

#Now fit an OU model that allows different sigma^2: 
out <- OUwie(tree,trait,model=c("OUMVA"),root.station=FALSE, plot.resid=FALSE)

# expecting mean 0, variance 1, with no differences between regimes 
print(out$Param.est)
