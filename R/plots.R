# File: plots.R
# Author: Carl Boettiger <cboettig@gmail.com>
# Date: 2011-11-16
# License: BSD


#' creates a color-coded plot of how well each tip fits it's predicted val
#' @param modelfit a fit from multiTypeOU (class multiOU)
#' @param n number of replicates to use in simulation
#' @param ... additional arguments to plot.phylo
#' @details 
#' this function is a proof of principle that could be made more generic
#' Note that this function also illustrates tree conversion with regime data.  
#' 
#' @export
tip_plot <- function(modelfit, n=100, ...){
  dat <- modelfit$data
  tree <- modelfit$tree

  reps <- sapply(1:n, function(i) t(simulate(modelfit)))

  m <- rowMeans(reps) 
  sder <- sapply(1:(dim(reps)[1]), function(i) sd(reps[i,]))
  colorcode <- dat[[1]]
  colorcode[ abs(dat - m) < 3*sder ] <- "orange"
  colorcode[ abs(dat - m) < 2*sder ] <- "green"
  colorcode[ abs(dat - m) < sder ] <- "blue"
  colorcode[ abs(dat - m) > 3*sder ] <- "red"
  tr <- convert(tree, colorcode)

  col <- tr$regimes[!is.na(tr$regimes)]
  plot(tr, ...)
  tiplabels(col=col, pch=19)
  list(tree=tr, tipcolors=col)
}








# These are data-format particular and could use updating.  
# Also consider ggplot-based utilities
# Meanwhile, these functions are not exported and plotting parameters is up to the user

subplot <- function(parname, par_dist, xlim=NULL, xlab=NULL, ...){
   colors <- c( rgb(0,1,0.5), rgb(0,0,1,.5), rgb(1,0,0,.5), rgb(1,1,0.5) )
   posterior <- vector("list", dim(par_dist)[2])
    id <- grep(parname, colnames(par_dist))
    for(i in id){
      posterior[[i]] <- density(par_dist[,i])
    }
    if(is.null(xlim[parname]))
      xlim <- c(min(sapply(posterior[id], function(P) P$x)),
                max(sapply(posterior[id], function(P) P$x)))
    else 
      xlim <- xlim[[parname]]

    if(is.null(xlab))
      xlab<-parse(text=parname)
    else
      xlab<-xlab[[parname]]

    plot(posterior[[id[1]]], xlab=xlab, xlim=xlim, ...)
    cindex <- 1
    for(i in id){
      polygon(posterior[[i]], col=colors[cindex])
      cindex <- cindex + 1
    }
}

plot_par <- function(par_dist, ...){
# Example:
#   plot_par(boots$null_par_dist, xlim=list(alpha=c(0,10)), cex=2)

  par_dist <- t(par_dist[!duplicated(par_dist),])
  par(mfrow=c(1,3))
  subplot("alpha", par_dist, main="Selection Strength", ...)
  subplot("sigma", par_dist, main="Trait Diversification", ...)
  subplot("theta", par_dist, main="Optimum", ...)
}



