# create a colorcoded plot of how well each tip is fitting it's predicted value
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


