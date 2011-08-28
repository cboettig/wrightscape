
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



