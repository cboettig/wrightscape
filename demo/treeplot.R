# CC0 

#To the extent possible under law, the author(s) have dedicated all copyright 
#and related and neighboring rights to this software to the public domain 
#worldwide. This software is distributed without any warranty. 

# For a copy of the CC0 Public Domain Dedication see,
# <http://creativecommons.org/publicdomain/zero/1.0/>.

# plot the phylogeny
rm(list=ls())
source("method2/method2_tools.R")

require(wrightscape)
require(geiger)
data(labrids)

row.names(ape.dat) <- ape.dat$Row.names

myphy <- paint_phy(ape.phy, ape.dat, list(c("Bolbometopon_muricatum", "Sparisoma_radians"), c("Chlorurus_sordidus", "Hipposcarus_longiceps")))


myphy$colors <- gsub("1", "blue",  myphy$colors)
myphy$colors <- gsub("2", "green",  myphy$colors)
myphy$colors <- gsub("3", "red",  myphy$colors)

cairo_pdf(file="phylo.pdf")
plot(myphy$phy, edge.color=myphy$colors, type="fan", show.tip.label=FALSE, edge.width=2)
dev.off()

