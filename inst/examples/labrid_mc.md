* Author: Carl Boettiger <cboettig@gmail.com>
* License: BSD 



```r
require(wrightscape)
require(snowfall)
require(ggplot2)
```







```r
data(labrids)
traits <- c("bodymass", "close", "open", "kt", "gape.y",  "prot.y", "AM.y", "SH.y", "LP.y")
```






```r
regimes <- two_shifts
```




Just a few processors, for debugging locally.


```r
sfInit(par=T, 4)    # for debugging locally
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




The main parallel loop fitting each model


```r
fits <- sfLapply(traits, function(trait){
	multi <- function(modelspec){ 
		out <- multiTypeOU(data = dat[[trait]], tree = tree, regimes = regimes, 
			    model_spec = modelspec, control = list(maxit=8000))
	      n <- length(levels(out$regimes))
	      Xo <- rep(out$Xo,n) 
	      loglik <- rep(out$loglik, n)
	      pars <- cbind(out$alpha, out$sigma, out$theta, Xo, loglik)
	      rownames(pars) <- levels(out$regimes)
	      colnames(pars) <- c("alpha", "sigma", "theta", "Xo", "loglik")
	      if(out$convergence != 0) # only return values if successful
      		pars[,] <- NA
	      pars
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
names(data) <- c("regimes", "param", "value", "model", "trait")
```




model likelihood



```r
ggplot(subset(data,  param=="loglik")) +
  geom_boxplot(aes(model, value)) +
  facet_wrap(~ trait, scales="free_y")
```

![plot of chunk unnamed-chunk-7](http://farm8.staticflickr.com/7114/6919293482_a5a6e34736_o.png) 


Parameter distributions of alpha parameter in model `alpha` (alphas vary) and `ou` (global).  



```r
ggplot(subset(data, param %in% c("alpha") 
       & model %in% c("alphas", "ou")),
       aes(model, value, fill=regimes)) +
  geom_bar(position="dodge") +  
  facet_wrap(~trait, scales="free_y")
```

![plot of chunk unnamed-chunk-8](http://farm6.staticflickr.com/5312/7065373031_0d9dec7c12_o.png) 


Parameter distribution of the sigma parameter in the brownie and bm models



```r
ggplot(subset(data, param %in% c("sigma") 
       & model %in% c("bm", "brownie")),
       aes(model, value, fill=regimes)) +
  geom_bar(position="dodge") +  
  facet_wrap(~trait, scales="free_y")
```

![plot of chunk unnamed-chunk-9](http://farm8.staticflickr.com/7277/6919294288_2686173e09_o.png) 



