rm(list=ls()) # clean workspace

script <- "method2_labrid_intra.R" 
#regime.names=c("wrasse", "pharyngeal", "intramandibular")
#regime.names=c("wrasse", "parrotfish")
regime.names=c("other", "intramandibular")

load("method2_labrid_intra.Rdat")



error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) != 
     length(lower) | length(lower) != length(upper))
  stop("vectors must be same length")
  arrows(x,y + upper, x, y - lower, 
         angle = 90, code = 3, length=length, ...)
}

alphas <- sapply(all_ouma, function(x){
  if(x$Diagnostic != "Arrived at a reliable solution")
    x$Param.est["alpha",] <- NA
  x$Param.est["alpha",]
})
colnames(alphas) <- X
alphas.se <- sapply(all_ouma, function(x){
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



sigma.sqs <- sapply(all_oumv, function(x){
  if(x$Diagnostic != "Arrived at a reliable solution")
    x$Param.est["sigma.sq",] <- NA
  x$Param.est["sigma.sq",]
})
colnames(sigma.sqs) <- X
sigma.sqs.se <- sapply(all_oumv, function(x){
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


require(socialR)
upload("test1.png test2.png", script=script, tags="phylogenetics")



