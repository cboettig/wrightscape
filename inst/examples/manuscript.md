* Author: Carl Boettiger <cboettig@gmail.com>
* License: BSD 






```r
require(wrightscape)
require(ggplot2)
require(reshape)
data(labrids)
```






```r
traits <- c("bodymass", "close", "open", "kt", "gape.y", "prot.y", 
    "AM.y")
regimes <- two_shifts
```







```r
nboot <- 100
cpu <- 16
```







```r
require(snowfall)
sfInit(parallel = TRUE, cpu = cpu)
```



```
R Version:  R version 2.15.0 (2012-03-30) 

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





Fit all models, then actually perform the model choice analysis for the chosen model pairs



```r
fits <- lapply(traits, function(trait) {
    multi <- function(modelspec) {
        multiTypeOU(data = dat[[trait]], tree = tree, regimes = regimes, model_spec = modelspec, 
            control = list(maxit = 8000))
    }
    bm2 <- multi(list(alpha = "fixed", sigma = "indep", theta = "global"))
    a2 <- multi(list(alpha = "indep", sigma = "global", theta = "global"))
    
    mc <- montecarlotest(bm2, a2, cpu = cpu, nboot = nboot)
})
```







```r
require(reshape2)
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
    # , Xo = df$Xo, converge=df$converge
    list(alpha = alpha, theta = theta, sigma = sigma)
}
dat <- melt(list(open = list(brownie = regroup(fits[[1]]$null_par_dist), 
    release = regroup(fits[[1]]$test_par_dist)), kt = list(brownie = regroup(fits[[2]]$null_par_dist), 
    release = regroup(fits[[2]]$test_par_dist))))
names(dat) <- c("regime", "value", "parameter", "model", "trait")
```








```r
ggplot(subset(dat, model == "release" & parameter == "alpha")) + 
    geom_boxplot(aes(model, value, fill = regime)) + facet_wrap(~trait, scale = "free_y")
```

![plot of chunk parameterplots](http://farm9.staticflickr.com/8147/7143478687_64dd96dba5_o.png) 

```r
ggplot(subset(dat, model == "brownie" & parameter == "sigma" & value < 
    2)) + geom_boxplot(aes(model, value, fill = regime)) + facet_wrap(~trait, 
    scale = "free_y")
```

![plot of chunk parameterplots](http://farm8.staticflickr.com/7070/6997392046_ac84be0a22_o.png) 


### Likelihood ratio distributions:



```r
lr_dat <- melt(list(open = list(release = fits[[1]]$test_dist, brownie = fits[[1]]$null_dist), 
    kt = list(release = fits[[2]]$test_dist, brownie = fits[[2]]$null_dist)))
open_LR <- 2 * (fits[[1]]$test$loglik - fits[[1]]$null$loglik)
kt_LR <- 2 * (fits[[2]]$test$loglik - fits[[2]]$null$loglik)
LR <- data.frame(value = c(open_LR, kt_LR), model = "ratio", trait = c("open", 
    "kt"))
names(lr_dat) <- c("value", "model", "trait")
```







```r
ggplot(subset(lr_dat, value > -10000 & value < 10000)) + geom_boxplot(aes(model, 
    value)) + facet_wrap(~trait, scale = "free_y") + geom_hline(data = LR, aes(yintercept = value))
```

![plot of chunk modelcomparisons](http://farm9.staticflickr.com/8026/6997392266_c9e253da17_o.png) 





```r
save(list = c("lr_dat", "dat", "fits"), file = "manuscript.rda")
```



