
data(geospiza)
attach(geospiza)


#---- PRINT RESULTS
fitContinuous(geospiza.tree, geospiza.data)
 
#---- STORE RESULTS 
brownFit <-  fitContinuous(geospiza.tree, geospiza.data)

aic.brown<-numeric(5)
for(i in 1:5) aic.brown[i]<-brownFit[[i]]$aic

#----------------------------------------------------
#   PHYLOGENETIC SIGNAL: FIT LAMBDA
#----------------------------------------------------

lambdaFit<-fitContinuous(geospiza.tree, geospiza.data, model="lambda")

# Compare likelihoods:

d.lambda<-numeric(5)
for(i in 1:5) d.lambda[i]=2*(lambdaFit[[i]]$lnl-brownFit[[i]]$lnl)

# Calculate p values assuming chi-squared distribution with 1 d.f.
p.lambda=pchisq(d.lambda, 1, lower.tail=FALSE)

aic.lambda<-numeric(5)
for(i in 1:5) aic.lambda[i]<-lambdaFit[[i]]$aic

#----------------------------------------------------
#    TIME PROPORTIONALITY: DELTA 
#---------------------------------------------------

deltaFit<-fitContinuous(geospiza.tree, geospiza.data, model="delta")

# Compare likelihoods:

d.delta<-numeric(5)
for(i in 1:5) d.delta[i]=2*(deltaFit[[i]]$lnl-brownFit[[i]]$lnl)

# Calculate p values assuming chi-squared distribution with 1 d.f.
p.delta=pchisq(d.delta, 1, lower.tail=FALSE)

aic.delta<-numeric(5)
for(i in 1:5) aic.delta[i]<-deltaFit[[i]]$aic

#----------------------------------------------------
#   SPECIATIONAL MODEL: KAPPA 
#---------------------------------------------------

kappaFit<-fitContinuous(geospiza.tree, geospiza.data, model="kappa")

# Compare likelihoods:

d.kappa<-numeric(5)
for(i in 1:5) d.kappa[i]=2*(kappaFit[[i]]$lnl-brownFit[[i]]$lnl)

# Calculate p values assuming chi-squared distribution with 1 d.f.
p.kappa=pchisq(d.kappa, 1, lower.tail=FALSE)

aic.kappa<-numeric(5)
for(i in 1:5) aic.kappa[i]<-kappaFit[[i]]$aic

#----------------------------------------------------
#   OU MODEL: ALPHA 
#---------------------------------------------------

ouFit<-fitContinuous(geospiza.tree, geospiza.data, model="OU")

# Compare likelihoods:

d.ou<-numeric(5)
for(i in 1:5) d.ou[i]=2*(ouFit[[i]]$lnl-brownFit[[i]]$lnl)

# Calculate p values assuming chi-squared distribution with 1 d.f.
p.ou=pchisq(d.ou, 1, lower.tail=FALSE)

aic.ou<-numeric(5)
for(i in 1:5) aic.ou[i]<-ouFit[[i]]$aic

#----------------------------------------------------
#   EARLY BURST MODEL: R 
#---------------------------------------------------

ebFit<-fitContinuous(geospiza.tree, geospiza.data, model="EB")

# Compare likelihoods:

d.eb<-numeric(5)
for(i in 1:5) d.eb[i]=2*(ebFit[[i]]$lnl-brownFit[[i]]$lnl)

# Calculate p values assuming chi-squared distribution with 1 d.f.
p.eb=pchisq(d.eb, 1, lower.tail=FALSE)

aic.eb<-numeric(5)
for(i in 1:5) aic.eb[i]<-ebFit[[i]]$aic

#----------------------------------------------------
#   COMPARE ALL MODELS
#---------------------------------------------------

# One way: use likelihood ratio test to compare all models to Brownian model

d.all<-cbind(d.lambda, d.delta, d.kappa, d.ou, d.eb)
p.all<-cbind(p.lambda, p.delta, p.kappa, p.ou, p.eb)

cat("Trait\tlambda\tdelta\tkappa\tou\teb\n")

for(i in 1:5) {
	cat("Tr", i, "\t");
	for(j in 1:5) {
		cat(round(d.all[i,j],2));
		if(p.all[i,j]<0.05) cat("*");
		if(p.all[i,j]<0.01) cat("*");
		if(p.all[i,j]<0.001) cat("*");
		cat("\t");
	}
	cat("\n");
}

# Another way: use AIC

aic.all<-cbind(aic.brown, aic.lambda, aic.delta, aic.kappa, aic.ou, aic.eb)
foo<-function(x) x-x[which(x==min(x))]
daic<-t(apply(aic.all, 1, foo))

rownames(daic)<-colnames(geospiza.data)
colnames(daic)<-c("Brownian", "Lambda", "Delta", "Kappa", "OU", "EB")

cat("Table of delta-aic values; zero - best model\n")
print(daic, digits=2)



