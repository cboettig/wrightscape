# knitr script
library(knitr)
pat_gfm()
render_gfm()
opts_knit$set(progress = TRUE, verbose = TRUE)
#opts_knit$set(upload.fun = socialR::flickr.url)
options(xtable.type = 'html')
opts_chunk$set(cache=TRUE, fig.path='figure/', dev='png', fig.width=8, fig.height=5, cache.path = 'cache/', par=TRUE)
opts_chunk$set(message=FALSE, warning=FALSE, comment=NA, tidy=FALSE)
options(device = function(width=5, height=5){ pdf(NULL, width=width, height=height) })

