rm(list=ls())
## Model Choice ##
load("labrid_power_d8bfe90.Rdat")
require(ggplot2)

dist_plot <- function(x, dat, compare, mode=c("density", "boxplot")){
  subdat <- subset(dat, comparison==compare & trait %in% x)
  if(mode == "density")
    ggplot(subset(subdat, type!="lr")) +
      geom_density(aes(value, fill=type),alpha=.8) +
      geom_vline(mapping = aes(yintercept = value),
                data = subset(subdat, type == "lr"), lty = 2) +
      facet_wrap( ~ trait, scale="free_y")
  else if (mode == "boxplot")
    ggplot(subset(subdat, type!="lr")) + 
      geom_boxplot(aes(type, value)) + 
     geom_hline(mapping = aes(yintercept = value),
                data = subset(subdat, type == "lr"), lty = 2) +
     facet_wrap( ~ trait, scale="free_y")
}

p1 <- dist_plot(c("kt", "open"), dat, "brownie_vs_alphas", "boxplot")
ggsave("model_choice.png", p1)


## parameters
rm(list=ls())
load("0a01a44.Rdat")
subdat <- subset(data, param %in% c("alpha") 
                 & model == "alphas" 
                 & trait %in% c("open", "kt"))
p2 <-  ggplot(subdat,aes(model, value, fill=regimes)) +
  geom_bar(position="dodge") +  
  facet_wrap(~trait, scales="free_y")
ggsave("parameters.png", p2)



## hmm 


load("f6cc2b1.Rdat")
require(Hmisc)
# Calculate the range for intellegent zooming in on summary stat values

subdat <- subset(data, param %in% c("alpha") 
                 & trait %in% c("kt", "close") 
                 & model %in% c("alphas"))
r <- cast(subdat, regimes ~ model ~ trait ~ param, smedian.hilow, conf.int=.5, na.rm=T)
upper <- sapply(c("alpha"), function(t) max(r[, , , t]))


p4 <-  ggplot(subdat, aes(model, value, fill=regimes)) + 
  stat_summary(fun.y=mean, geom="bar", position="dodge", alpha=.5) + # add bars for some extra ink...
  stat_summary(fun.data=median_hilow, geom="pointrange", aes(color=regimes), 
  position = position_dodge(width=0.90), conf.int=.5) +
  facet_grid(param ~ trait, scales = "free_y") + 
  coord_cartesian(ylim=c(0,upper["alpha"]), wise=TRUE) +
  opts(title = "alpha")


