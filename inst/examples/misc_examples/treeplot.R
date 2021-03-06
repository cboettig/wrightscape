# CC0 

#To the extent possible under law, the author(s) have dedicated all copyright 
#and related and neighboring rights to this software to the public domain 
#worldwide. This software is distributed without any warranty. 

# For a copy of the CC0 Public Domain Dedication see,
# <http://creativecommons.org/publicdomain/zero/1.0/>.

# plot the phylogeny
rm(list=ls())
source("method2_tools.R")

myphy <- paint_phy(ape.phy, ape.dat, list(c("Bolbometopon_muricatum", "Sparisoma_radians"), c("Chlorurus_sordidus", "Hipposcarus_longiceps")))

cairo_pdf(file="phylo.pdf")
plot(myphy$phy, edge.color=myphy$colors, type="fan", show.tip.label=FALSE, edge.width=2)
dev.off()

