# plot the phylogeny
rm(list=ls())
source("method2_tools.R")

myphy <- paint_phy(ape.phy, ape.dat, list(c("Bolbometopon_muricatum", "Sparisoma_radians"), c("Chlorurus_sordidus", "Hipposcarus_longiceps")))

cairo_pdf(file="phylo.pdf")
plot(myphy$phy, edge.color=myphy$colors, type="fan", show.tip.label=FALSE, edge.width=2)
dev.off()

