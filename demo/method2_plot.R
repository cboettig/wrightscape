rm(list=ls()) # clean workspace
load("method2_parrotfish.Rdat")
#load("method2_labrid.Rdat")



error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) != 
     length(lower) | length(lower) != length(upper))
  stop("vectors must be same length")
  arrows(x,y + upper, x, y - lower, 
         angle = 90, code = 3, length=length, ...)
}


regime.names=c("other", "intramandibular")

alphas <- sapply(all_oumva, function(x){
  if(x$Diagnostic != "Arrived at a reliable solution")
    x$Param.est["alpha",] <- NA
  x$Param.est["alpha",]
})
colnames(alphas) <- X
alphas.se <- sapply(all_oumva, function(x){
  if(x$Diagnostic != "Arrived at a reliable solution")
    x$Param.SE["alpha",] <- NA
  x$Param.SE["alpha",]
})

#### Plot alphas ###
png("test1.png", width=600)
  bars <- barplot(alphas, beside=T, main="alphas", 
  legend.text=regime.names,
  ylim=c(0, max(alphas+alphas.se, na.rm=T)))
  error.bar(bars, alphas, alphas.se)
dev.off()



sigma.sqs <- sapply(all_oumva, function(x){
  if(x$Diagnostic != "Arrived at a reliable solution")
    x$Param.est["sigma.sq",] <- NA
  x$Param.est["sigma.sq",]
})
colnames(sigma.sqs) <- X
sigma.sqs.se <- sapply(all_oumva, function(x){
  if(x$Diagnostic != "Arrived at a reliable solution")
    x$Param.SE["sigma.sq",] <- NA
  x$Param.SE["sigma.sq",]
})
#### Plot sigma.sq ###
png("test2.png", width=600)
  bars <- barplot(sigma.sqs, beside=T,
                  main=paste(expression(sigma^2)),
                  legend.text=regime.names,
                 ylim=c(0, max(sigma.sqs+sigma.sqs.se, na.rm=T)))
  error.bar(bars, sigma.sqs, sigma.sqs.se)
dev.off()

thetas <- sapply(all_oumva, function(x){
  if(x$Diagnostic != "Arrived at a reliable solution")
    x$theta[,"Estimate"] <- NA
  x$theta[,"Estimate"]
})
colnames(thetas) <- X
thetas.se <- sapply(all_oumva, function(x){
  if(x$Diagnostic != "Arrived at a reliable solution")
    x$theta[,"SE"] <- NA
  x$theta[,"SE"]
})
#### Plot thetas ###
png("test3.png", width=600)
  bars <- barplot(thetas, beside=T, main="thetas", 
                 legend.text=regime.names, 
                 ylim=c(min(thetas-thetas.se, na.rm=T),
                 max(thetas+thetas.se, na.rm=T)))
  error.bar(bars, thetas, thetas.se)
  legend
dev.off()


require(socialR)
upload("test*.png", script="method2_parrotfish.R", tags="phylogenetics")

load("method2_labrid.Rdat")

png("phylo.png", 600, 600)
input <- paint_phy(ape$phy, traits,  
         list(c("Bolbometopon_muricatum", "Sparisoma_radians"), 
              c("Chlorurus_sordidus", "Hipposcarus_longiceps")))
#input <- paint_phy(ape$phy, traits,  c("Chlorurus_sordidus", "Hipposcarus_longiceps"))
dev.off()
require(socialR)
flickr("phylo.png", tag="phylogenetics")

