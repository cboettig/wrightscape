* File: labrid.R
* Author: Carl Boettiger <cboettig@gmail.com>
* License: BSD 





```r
require(wrightscape)
require(snowfall)
require(ggplot2)
require(reshape2)
require(Hmisc)
```








```r
data(labrids)
traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")
regimes <- two_shifts
```




Setup parallel environment



```r
sfInit(par=T, 4)   
```



```
R Version:  R version 2.14.1 (2011-12-22) 

```



```r
sfLibrary(wrightscape)
```



```
Library wrightscape loaded.
```



```r
sfExportAll()
```




Bootstrap the fits by replication (bootstraps parameter values without bootstrapping the likelihood ratios)




```r
fits <- sfLapply(traits, function(trait){
  multi <- function(modelspec, reps = 100){
    m <- multiTypeOU(data = dat[trait], tree = tree, regimes = regimes, 
  		     model_spec = modelspec, 
		     control = list(maxit=8000)
		    ) 
    replicate(reps, bootstrap(m))
  }
  bm <- multi(list(alpha = "fixed", sigma = "global", theta = "global")) 
  ou <- multi(list(alpha = "global", sigma = "global", theta = "global")) 
  bm2 <- multi(list(alpha = "fixed", sigma = "indep", theta = "global")) 
  a2  <- multi(list(alpha = "indep", sigma = "global", theta = "global")) 
  t2  <- multi(list(alpha = "global", sigma = "global", theta = "indep")) 
	list(bm=bm,brownie=bm2, ou=ou, ouch=t2, alphas=a2)
})
```




Reformat and label data for plotting



```r
names(fits) <- traits  # each fit is a different trait (so use it for a label)
data <- melt(fits)
names(data) <- c("regimes", "param", "rep", "value", "model", "trait")
```







```r
subdat <- subset(data, param %in% c("alpha") 
                 & model %in% c("alphas") 
                 & value < 20)
r <- cast(subdat, regimes ~ model ~ trait ~ param, smedian.hilow, conf.int=.5, na.rm=T)
upper <- sapply(c("alpha"), function(t) max(r[, , , t]))
```





Generate a plot of the relevant parameters



```r
p0 <- 
ggplot(subdat, aes(model, value, fill=regimes)) + 
  stat_summary(fun.data=median_hilow, geom="pointrange", aes(color=regimes), 
  position = position_dodge(width=0.90), conf.int=.5) +
  facet_grid(param ~ trait, scales = "free_y") + 
  coord_cartesian(ylim=c(0,upper["alpha"]), wise=TRUE) +
  opts(title = "alpha")
p0
```

![plot of chunk unnamed-chunk-7](http://farm8.staticflickr.com/7256/6922704638_684d116bb8_o.png) 






```r
save(list=ls, file = "~/public_html/data/labrids_par_bootstrap.rda")
```



```
Error: invalid first argument
```






Additional plots, just for variety



```r
range <- cast(data, regimes ~ model ~ trait ~ param, smedian.hilow,
              conf.int=.5, na.rm=T)
upper <- sapply(c("sigma", "alpha"), function(t) max(range[, , , t]))
p1 <- ggplot(subset(data,  param=="loglik")) + 
      geom_boxplot(aes(model, value)) +
      facet_wrap(~ trait, scales="free_y")
p1
```

![plot of chunk unnamed-chunk-9](http://farm8.staticflickr.com/7253/7066357345_5b7b949d07_o.png) 

```r
p2 <- ggplot(subset(data, param %in% c("sigma", "alpha") & model != "bm"),
              aes(model, value, fill=regimes)) + 
       stat_summary(fun.y=mean, geom="bar", position="dodge") + 
       stat_summary(fun.data=median_hilow, geom="pointrange",
       aes(color=regimes), 
       position = position_dodge(width=0.90), conf.int=.5) +
       scale_y_log() + 
       facet_grid(param ~ trait, scales = "free_y") + 
       coord_cartesian(ylim=c(0,max(upper)), wise=TRUE)
p2
```

![plot of chunk unnamed-chunk-9](http://farm6.staticflickr.com/5333/7066357573_3e71012e4e_o.png) 

```r
p3 <- ggplot(subset(data, param %in% c("sigma") ), 
              aes(model, value, fill=regimes)) + 
       stat_summary(fun.y=mean, geom="bar", position="dodge", alpha=.5) + 
       stat_summary(fun.data=median_hilow, geom="pointrange",
       aes(color=regimes), 
       position = position_dodge(width=0.90), conf.int=.5) +
       facet_grid(param ~ trait, scales = "free_y") + 
       coord_cartesian(ylim=c(0,upper["sigma"]), wise=TRUE) +  
       opts(title="sigma")
p3
```

![plot of chunk unnamed-chunk-9](http://farm8.staticflickr.com/7251/6920279016_163f5c778f_o.png) 

```r
p4 <- ggplot(subset(data, param %in% c("alpha") ), 
              aes(model, value, fill=regimes)) + 
       stat_summary(fun.y=mean, geom="bar", position="dodge", alpha=.5) +
       stat_summary(fun.data=median_hilow, geom="pointrange",
       aes(color=regimes), 
       position = position_dodge(width=0.90), conf.int=.5) +
       facet_grid(param ~ trait, scales = "free_y") + 
       coord_cartesian(ylim=c(0,upper["alpha"]), wise=TRUE) +
        opts(title = "alpha")
p4
```

![plot of chunk unnamed-chunk-9](http://farm6.staticflickr.com/5450/7066358019_01ae5deb6a_o.png) 




```r
save(list=ls, file = "~/public_html/data/labrids_par_bootstrap.rda")
```



```
Error: invalid first argument
```






