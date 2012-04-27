* Author: Carl Boettiger <cboettig@gmail.com>
* License: BSD 




```r
require(wrightscape)
```



```
## Loading required package: wrightscape
```



```
## Loading required package: mcmcTools
```



```
## Loading required package: ape
```



```
## 
## Attaching package: 'wrightscape'
## 
```



```
## The following object(s) are masked from 'package:mcmcTools':
## 
##     getParameters
## 
```



```r
require(ggplot2)
```



```
## Loading required package: ggplot2
```



```
## Loading required package: reshape
```



```
## Loading required package: plyr
```



```
## 
## Attaching package: 'reshape'
## 
```



```
## The following object(s) are masked from 'package:plyr':
## 
##     rename, round_any
## 
```



```
## Loading required package: grid
```



```
## Loading required package: proto
```



```r
require(reshape)
data(labrids)
```






```r
traits <- c("open", "kt")
regimes <- two_shifts
```







```r
nboot <- 100
cpu <- 16
```







```r
require(snowfall)
```



```
## Loading required package: snowfall
```



```
## Loading required package: snow
```



```r
sfInit(parallel = TRUE, cpu = cpu)
```



```
## R Version:  R version 2.14.1 (2011-12-22) 
## 
```



```
## snowfall 1.84 initialized (using snow 0.3-8): parallel execution on 16 CPUs.
## 
```



```r
sfLibrary(wrightscape)
```



```
## Library wrightscape loaded.
```



```
## Library wrightscape loaded in cluster.
## 
```



```
## Warning message: 'keep.source' is deprecated and will be ignored
```



```r
sfExportAll()
```





Fit all models, then actually perform the model choice analysis for the chosen model pairs



```r
fits <- lapply(traits, function(trait) {
    multi <- function(modelspec) {
        multiTypeOU(data = dat[[trait]], tree = tree, regimes = regimes, 
            model_spec = modelspec, control = list(maxit = 8000))
    }
    bm2 <- multi(list(alpha = "fixed", sigma = "indep", theta = "global"))
    a2 <- multi(list(alpha = "indep", sigma = "global", theta = "global"))
    
    mc <- montecarlotest(bm2, a2, cpu = cpu, nboot = nboot)
})
```





Clean up the data



```r
require(reshape2)
```



```
## Loading required package: reshape2
```



```
## 
## Attaching package: 'reshape2'
## 
```



```
## The following object(s) are masked from 'package:reshape':
## 
##     colsplit, melt, recast
## 
```



```r
require(ggplot2)
```




### Parameter distributions



```r
regroup <- function(df) {
    df <- as.data.frame(t(df))
    alpha <- df[c("alpha1", "alpha2", "alpha3")]
    sigma <- df[c("sigma1", "sigma2", "sigma3")]
    theta <- df[c("theta1", "theta2", "theta3")]
    names(alpha) <- levels(fits[[1]]$null$regimes)
    names(sigma) <- levels(fits[[1]]$null$regimes)
    names(theta) <- levels(fits[[1]]$null$regimes)
    #, Xo = df$Xo, converge=df$converge
    list(alpha = alpha, theta = theta, sigma = sigma)
}
dat <- melt(list(open = list(brownie = regroup(fits[[1]]$null_par_dist), 
    release = regroup(fits[[1]]$test_par_dist)), kt = list(brownie = regroup(fits[[2]]$null_par_dist), 
    release = regroup(fits[[2]]$test_par_dist))))
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```
## Using  as id variables
```



```r
names(dat) <- c("regime", "value", "parameter", 
    "model", "trait")
```








```r
ggplot(subset(dat, model == "release" & parameter == 
    "alpha")) + geom_boxplot(aes(model, value, fill = regime)) + 
    facet_wrap(~trait, scale = "free_y")
```

![plot of chunk unnamed-chunk-8](http://farm9.staticflickr.com/8154/6971407232_7e56beaa36_o.png) 

```r
ggplot(subset(dat, model == "brownie" & parameter == 
    "sigma" & value < 2)) + geom_boxplot(aes(model, value, 
    fill = regime)) + facet_wrap(~trait, scale = "free_y")
```

![plot of chunk unnamed-chunk-8](http://farm8.staticflickr.com/7110/7117484781_ff3950f264_o.png) 


### Likelihood ratio distributions:



```r
lr_dat <- melt(list(open = list(release = fits[[1]]$test_dist, 
    brownie = fits[[1]]$null_dist), kt = list(release = fits[[2]]$test_dist, 
    brownie = fits[[2]]$null_dist)))
open_LR <- 2 * (fits[[1]]$test$loglik - fits[[1]]$null$loglik)
kt_LR <- 2 * (fits[[2]]$test$loglik - fits[[2]]$null$loglik)
LR <- data.frame(value = c(open_LR, kt_LR), model = "ratio", 
    trait = c("open", "kt"))
names(lr_dat) <- c("value", "model", "trait")
```







```r
ggplot(subset(lr_dat, value > -10000 & value < 
    10000)) + geom_density(aes(value, fill = model), alpha = 0.5) + 
    facet_wrap(~trait, scale = "free_x") + geom_vline(data = LR, 
    aes(xintercept = value))
```

![plot of chunk unnamed-chunk-10](http://farm9.staticflickr.com/8022/7117485015_969c66685f_o.png) 




```r
save(list = c("lr_dat", "dat", "fits"), file = "manuscript.rda")
```



