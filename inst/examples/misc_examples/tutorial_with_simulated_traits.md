# A wrightscape tutorial using simulated traits



Load package, along with parallelization, plotting, and data manipulation packages. Then load the data.  




```r
require(wrightscape)
require(snowfall)
require(ggplot2)
require(reshape)
data(labrids)
```





rename the regimes less technically



```r
levels(pharyngeal) = c("wrasses", "parrotfish")
regime <- pharyngeal
```




Create a dataset by simulation, where the parrotfish all have lower alpha



```r
a1_spec  <- list(alpha = "indep", sigma = "global", theta = "global")
a1 <- multiTypeOU(data = dat["close"], tree = tree, regimes = regime, 
	     model_spec = a1_spec,  control = list(maxit=5000))
```




assign expected "startree" standard deviations: wrasses 1.25, parrotfish 5



```r
a1$alpha[1] <- 10
a1$alpha[2] <- 1e-10
a1$sigma <- c(5, 5)  
a1$theta <- c(0,0)   
dat[["constraint release"]] <-simulate(a1)[[1]]
```




Check out the variance in the relative groups -- it should be larger in parrotfish
in an extreme example, but need not be in principle, due to the phylogeny. 



```r
testcase <- dat[["constraint release"]]
testcase[regime =="wrasses" & !is.na(testcase) & testcase != 0] -> lowvar
testcase[regime !="wrasses" & !is.na(testcase) & testcase != 0] -> highvar
print(c(var(lowvar), var(highvar)))
```



```
[1]  1.471 10.988
```




We can repeat the whole thing with a model based on differnt sigmas, to make sure 


```r
s1_spec  <- list(alpha = "global", sigma = "indep", theta = "global")
s1 <- multiTypeOU(data = dat["close"], tree = tree, regimes = pharyngeal, 
	     model_spec = s1_spec,  control = list(maxit=5000))
names(s1$sigma) <- levels(regime)
s1$sigma[1] <- sqrt(2*5*1.25)
s1$sigma[2] <- sqrt(2*5*5)
s1$alpha <- c(5, 5)  # We can keep those parameters estimated from data or update them
a1$theta <- c(0,0)   
dat[["faster evolution"]] <-simulate(s1)[[1]]
testcase <- dat[["faster evolution"]]
testcase[regime == "wrasses" & !is.na(testcase) & testcase != 0 ] -> lowvar
testcase[regime != "wrasses" & !is.na(testcase) & testcase != 0 ] -> highvar
print(c(var(lowvar), var(highvar)))
```



```
[1] 1.393 5.884
```




Now we have a trait where change in alpha is responsible, 
and one in which sigma change is responsible. 
Can we correctly identify each??



```r
traits <- c("constraint release", "faster evolution")
```






```r
fits <- lapply(traits, function(trait){
  # declare function for shorthand
  multi <- function(modelspec, reps = 20){
    m <- multiTypeOU(data = dat[[trait]], tree = tree, regimes = pharyngeal, 
  		     model_spec = modelspec, 
		     control = list(maxit=5000)
		    ) 
    replicate(reps, bootstrap(m))
  }
  bm <- multi(list(alpha = "fixed", sigma = "indep", theta = "global"))
  s1 <- multi(list(alpha = "global", sigma = "indep", theta = "global")) 
  a1  <- multi(list(alpha = "indep", sigma = "global", theta = "global")) 
  s2 <- multi(list(alpha = "global", sigma = "indep", theta = "indep")) 
  a2  <- multi(list(alpha = "indep", sigma = "global", theta = "indep")) 
  list(bm=bm, s1=s1, a1=a1, s2=s2, a2=a2)
})
```




Reformat and label data for plotting



```r
names(fits) <- traits  # each fit is a different trait (so use it for a label)
data <- melt(fits)
names(data) <- c("regimes", "param", "rep", "value", "model", "trait")
```






```r
save(list=ls(), file="tutorial_with_simulated_traits.rda")
```





model likelihood


```r
p1 <- ggplot(subset(data,  param=="loglik")) + 
      geom_boxplot(aes(model, value)) +
      facet_wrap(~ trait, scales="free_y")
p1
```

![plot of chunk unnamed-chunk-11](http://farm8.staticflickr.com/7065/7068904203_ac4675cc7a_o.png) 




```r
p2 <-  ggplot(subset(data, param %in% c("sigma", "alpha")), aes(model, value, fill=regimes)) + 
       stat_summary(fun.data=mean_sdl, geom="pointrange", aes(color=regimes), 
		    position = position_dodge(width=0.90)) +
       scale_y_log() + 
       facet_grid(param ~ trait, scales = "free_y")
```



```
Error: could not find function "scale_y_log"
```



```r
p2
```



```
Error: object 'p2' not found
```



```r
p3 <- ggplot(subset(data, param %in% c("sigma") )) +
      geom_boxplot(aes(model, value, fill=regimes)) + 
      facet_wrap(trait ~ param, scales = "free_y") 
p3
```

![plot of chunk unnamed-chunk-12](http://farm6.staticflickr.com/5448/6922823962_e9f978c25e_o.png) 

```r
p4 <- ggplot(subset(data, param %in% c("alpha")  )) +
      geom_boxplot(aes(model, value, fill=regimes)) + 
      facet_wrap(trait ~ param, scales = "free_y") 
p4
```

![plot of chunk unnamed-chunk-12](http://farm6.staticflickr.com/5446/7068904775_d31fcda8b3_o.png) 




```r
save(list=ls(), file="tutorial_with_simulated_traits.rda")
```





