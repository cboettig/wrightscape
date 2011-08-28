# evolution plots

require(wrightscape)
require(pmc)
load("5835164312.Rdat")
names(sigmas$data) # note sigmas is also in boots$null


#c("close", "open", "gape.y", "prot.y")
sapply(labrid$data[c(1:14)], function(trait){
  var(trait[sigmas$regimes=="intramandibular"], na.rm=T)/var(trait[sigmas$regimes=="other"], na.rm=T) 
  })


sapply(labrid$data[c(3,7:14)], function(trait){
  var(trait[sigmas$regimes=="intramandibular"], na.rm=T)/var(trait[sigmas$regimes=="other"], na.rm=T) 
  })

sapply(labrid$data[c(3,7:14)], function(trait){
  mean(trait[sigmas$regimes=="intramandibular"], na.rm=T)/mean(trait[sigmas$regimes=="other"], na.rm=T) 
  })



info_plot <- function(null, test, nullname="Null Model", testname="Test Model"){
  options(digits=2)
  plot(0,0, type="n", xaxt="n", yaxt="n", xlab="", ylab="")

  
  tmp <- getParameters(null)
  tmp <- tmp[!duplicated(tmp)]

  text(-.2, .7, nullname, cex=1.5)
  text(-1, .5, paste(names(tmp), collapse="\t"), pos=4)
  text(-1, .4, paste(formatC(tmp, format="fg", digits=4),
       collapse="\t"), pos=4) 

  tmp <- getParameters(test)
  tmp <- tmp[!duplicated(tmp)]


  text(-.2, 0, testname, cex=1.5)
  text(-1, -.2, paste(names(tmp), collapse="\t"), pos=4) 
  text(-1, -.3, paste(formatC(tmp, format="fg", digits=4),
      collapse="\t"), pos=4) 

  text(-1, -.6, paste("Likelihood Ratio = ", 
       prettyNum(-2*(loglik(null)-loglik(test)))),
       cex=1.5, pos=4) 

}
info_plot(sigmas, alphas, "Disparity-driven", "Constraint-driven")

source("../../pmc/R/treeformats.R")
ape <- convert(sigmas$tree, intramandibular)


png("parrotfish_tree.png", height=1000, width=1000)
plot(ape, edge.color=treepalette(ape, custom=c("green", "blue")), edge.width=5,cex=1.5)
legend("bottomleft", c("Intramandibular Innovation", "Others"), pch=15, col=c("green", "blue"), cex=3)
dev.off()


#plot_par(boots$null_par_dist)
#plot_par(boots$test_par_dist, xlim=list(alpha=c(0,6), sigma=c(0,.6), theta=c(-0.5,2)))




