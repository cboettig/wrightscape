# plot all the power distributions as 95\% confidence intervals

rm(list=ls())
#load("labrid_power_6befae4.Rdat")
load("labrid_power_4f97994.Rdat")
names(fits) <- traits
dat <- melt(fits)

names(dat) <- c("value", "type", "trait")
dat[dat$type == 1, 2] <- "null"
dat[dat$type == 2, 2] <- "test"
dat[dat$type == 3, 2] <- "lr"


require(ggplot2)

p1 <- ggplot(subset(dat, abs(value) < 1e3)) + 
      geom_boxplot(aes(type, value)) +
      facet_wrap(~ trait, scales="free_y")
print(p1)
