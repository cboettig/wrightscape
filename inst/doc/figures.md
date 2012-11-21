





```r
px
```

![plot of chunk alpha_v_sigma](figure/alpha_v_sigma1.pdf) 

```r
py
```

![plot of chunk alpha_v_sigma](figure/alpha_v_sigma2.pdf) 

```r
pw
```

![plot of chunk alpha_v_sigma](figure/alpha_v_sigma3.pdf) 

```r
pz
```

![plot of chunk alpha_v_sigma](figure/alpha_v_sigma4.pdf) 





```r
load("wrightscape.rda")
```








```r
ggplot(subset(lr_dat, value > -10000 & value < 10000)) + geom_density(aes(value, 
    color = model), alpha = 0.8, color = c("blue", "red")) + facet_wrap(~trait, 
    scale = "free_y") + geom_vline(data = LR, aes(xintercept = value)) + ylab("likelihood ratio")
```

```
## Error: When _setting_ aesthetics, they may only take one value. Problems:
## colour
```





