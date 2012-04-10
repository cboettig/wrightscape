library(knitr)
pat_gfm()
render_gfm()
opts_knit$set(progress = TRUE, verbose = TRUE)
opts_knit$set(upload.fun = socialR::flickr.url)
options(xtable.type = 'html')
opts_chunk$set(cache=TRUE, fig.path='figure/', dev='Cairo_png',
  fig.width=6, fig.height=5, cache.path = 'cache-local/', par=TRUE)
opts_chunk$set(message=FALSE, warning=FALSE, comment=NA, tidy=FALSE)

