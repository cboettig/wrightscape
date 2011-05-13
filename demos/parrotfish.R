# parrotfish.R
require(wrightscape)
require(pmc)
require(socialR)
tag="phylogenetics wrightscape labrids"

source("parrotfish_data.R")
source("loop_models_traits_regimes.R")

model_list <- list("brown", "hansen", "ouch", "brownie", "wright", "release_constraint")
regime_list <-  list(intramandibular=intramandibular)

#test <- fit_all(model_list, labrid$data, regime_list, labrid$tree)

test <- fit_all(model_list[4:6], labrid$data[10:13], regime_list, labrid$tree)

conv(test)

#test$fits[[1]][[1]][[3]]$optim_output$par


## Summary matrices
likmat <- llik_matrix(test)
lliks <- sort_lik(likmat)
release_alphas <- alpha_traits(test) 
wright_alphas <- alpha_traits(test, model_name="wright") 


#### Plots

# Body mass, ratios 
png("parrotfish_ratios.png", width=480*1.5, height=480*1.5)
par(mfrow=c(2,2))
for(i in c(3, 7:9)){
  barplot(lliks[[i]], horiz=T, main=names(lliks)[i])
}
dev.off()
# Body mass, corrected masses
png("parrotfish_mass_y.png", width=480*1.5, height=480*1.5)
par(mfrow=c(2,2))
for(i in c(3, 12:14)){
  barplot(lliks[[i]], horiz=T, main=names(lliks)[i])
}
dev.off()
# Body mass and uncorrected masses
png("parrotfish_mass_x.png", width=480*1.5, height=480*1.5)
par(mfrow=c(2,2))
for(i in c(3:6)){
  barplot(lliks[[i]], horiz=T, main=names(lliks)[i])
}
dev.off()
# Corrected and uncorrected lengths
png("parrotfish_lengths.png", width=480*1.5, height=480*1.5)
par(mfrow=c(2,2))
for(i in c(1,2,10,11)){
  barplot(lliks[[i]], horiz=T, main=names(lliks)[i])
}
dev.off()


png("release_alphas.png", width=480*2)
 barplot(log(release_alphas[,1]/release_alphas[,2]))
dev.off()


png("wright_alphas.png", width=480*2)
 barplot(log(wright_alphas[,1]/wright_alphas[,2]))
dev.off()

flickr(files="parrotfish*.png", tag=tag)
flickr(files="*alphas.png", tag=tag)








## sometimes hansen returns a silly large value of alpha, which will break the generalized likelihood function (too stiff)
 #h <- fit("hansen", labrid$data[13], intramandibular, labrid$tree, .01, .01)
 # w <- fit("release_constraint", labrid$data[13], intramandibular, labrid$tree, h@sqrt.alpha^2, h@sigma)



