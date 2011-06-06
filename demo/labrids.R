# labrid example

require(wrightscape)
require(pmc)
require(socialR)
tag="phylogenetics wrightscape labrids"

# This data has not been released
path = "../data/labrids/"
labrid_tree <- read.nexus(paste(path, "labrid_tree.nex", sep=""))
#fin_data <-read.csv(paste(path,"labrid.csv", sep=""))
diet_data <- read.csv(paste(path,"labriddata_parrotfish.csv", sep=""))

## Unfortunately this "phylogenetically accurate" size correction assumes 
source("correct_labrid_data.R")
source("paint_labrid_regimes.R")
source("loop_models_traits_regimes.R")

model_list <- list("brown", "hansen", "ouch", "brownie", "wright", "release_constraint")
regime_list <-  list(pharyngeal=pharyngeal, intramandibular=intramandibular, two_shifts=two_shifts)


test <- results(model_list[c(1,2,4,6)], labrid$data[10:13], regime_list, labrid$tree)

mm <- llik_matrix(model_list[c(1,2,4,6)], names(regime_list), j=1)
barplot(mm, col=c("thistle", "khaki", "palegreen"), horiz=TRUE, beside=TRUE, legend=TRUE)

#labrid_results <- results(model_list, labrid$data[1:3], regime_list, labrid$tree)



#argument is: (traits,regime,model) i.e. (j,k,i)

llik_matrix <- function(models, regimes, j, results){
# 
# models: a list of models (character vector)
# regimes: a list of regimes (character vector)

  M <- matrix(NA, nrow=length(regimes), ncol=length(models)) 
  for(i in 1:length(models)){
    for(k in 1:length(regimes)){
      a <- results[[j]][[k]][[i]]
      if (is(a, "try-error")){
        M[k,i] <- NA
      } else if(is(a, "hasentree") | is(a, "browntree")){
        M[k,i] <- a@loglik
      } else {
        M[k,i] <- loglik(a)
      }
    }
    rownames(M) <- regimes
    colnames(M) <- models
  }
  M
}



for(i in 1:length(labrid_results)){
mm <- llik_matrix(model_list, names(regime_list), j=i)
png(paste(i,".png"))
barplot(mm, col=c("thistle", "khaki", "palegreen"), horiz=TRUE, 
        beside=TRUE, legend=TRUE)
dev.off()
}

png("gape.png")
barplot(gape, col=c("thistle", "khaki", "palegreen"), horiz=TRUE, 
        beside=TRUE, main="gape", legend=TRUE, xlim=c(-90,-65))
dev.off()



alpha_matrix <- function(model, regime_list, regime_k, trait_i){
  i <- trait_i
  k <- regime_k
  regime <- regime_list[[k]]
  nregimes <- length(levels(regime))
  M <- matrix(NA, nrow=nregimes, ncol=length(model)) 
  for(i in 1:length(model)){
      a <- labrid_results[[j]][[k]][[i]]
    for(k in 1:nregimes){
      if (is(a, "try-error") | is(a, "browntree")){
        print(class(a))
        M[k,i] <- NA
      } else if(is(a, "hansentree")){
        M[k,i] <- (a@sqrt.alpha)^2
      } else if(is(a, "multiOU")) {
        M[k,i] <- a$alpha[k]
      }
    }
    colnames(M) <- model
    rownames(M) <- levels(regime)
  }
  M
}
j <- 1; k <- 3

for(j in 1:3){
gape_alpha <- alpha_matrix(model_list[c(1,2,6)], regime_list=regime_list, regime_k=3, trait=j)

png(paste(j,"_alpha.png", sep=""))
barplot(gape_alpha, col=c("thistle", "khaki", "palegreen"), horiz=TRUE, 
        beside=TRUE, main="gape", legend=TRUE)
dev.off()
}









dummy <- function(i){
	trait_name <- names(labrid$data)[i]	
	trait <- labrid$data[i]


	ou1 <- hansen(trait, labrid$tree, regime=labrid$noregimes, .01, .01)
	ou2_phar <- hansen(trait, labrid$tree, regime=labrid$regimes, .01, .01)
  ou2_intra <- hansen(trait, labrid$tree, regime=intramandibular, .01, .01)
  ou3 <- hansen(trait, labrid$tree, regime=two_shifts, ou2_phar@sqrt.alpha, sigma=ou2_phar@sigma )

  ## Note that this will converge poorly with the .01, .01 starting conditions
  #  ou3 <- hansen(trait, labrid$tree, regime=two_shifts, .01, .01 )
  #  loglik(ou3) - loglik(ou2_phar)



  ## Compare intramandibular and pharyngeal joint paintings
  ##  Note that we use the more restricted estimates from ouch to seed the search

#  pharyngeal vs intramadibular regimes
  ouch_phar <- ouch(trait, labrid$tree, regime=labrid$regimes, 
                    alpha=(ou2_phar@ sqrt.alpha)^2, sigma=ou2_phar@sigma)
	brownie_phar <- brownie(trait, labrid$tree, regime=labrid$regimes,
                          sigma=ou2_phar@sigma)
	wright_phar <- wright(trait, labrid$tree, regime=labrid$regimes,
                        alpha=(ou2_phar@sqrt.alpha)^2, sigma=ou2_phar@sigma)
  release_phar <- release_constraint(trait, labrid$tree, regime=labrid$regimes,
                                     alpha=(ou2_phar@sqrt.alpha)^2,
                                     sigma=ou2_phar@sigma)

  ouch_intra <- ouch(trait, labrid$tree, regime=intramandibular,
                     alpha=(ou2_phar@sqrt.alpha)^2, sigma=ou2_phar@sigma)
	brownie_intra <- brownie(trait, labrid$tree, regime=intramandibular,
                           sigma=ou2_intra@sigma)
	wright_intra <- wright(trait, labrid$tree, regime=intramandibular,
                         alpha=(ou2_intra@sqrt.alpha)^2, sigma=ou2_intra@sigma)
  release_intra <- release_constraint(trait, labrid$tree,
                                      regime=intramandibular, 
                                      alpha=(ou2_intra@sqrt.alpha)^2, 
                                      sigma=ou2_intra@sigma)

  ouch_twoshifts <- ouch(trait, labrid$tree, regime=two_shifts,
                        alpha=(ou3@sqrt.alpha)^2, sigma=ou3@sigma)
	brownie_twoshifts <- brownie(trait, labrid$tree, regime=two_shifts,
                               sigma=ou3@sigma)
	wright_twoshifts <- wright(trait, labrid$tree, regime=two_shifts,
                              alpha=(ou3@sqrt.alpha)^2, sigma=ou3@sigma)
  release_twoshifts <- release_constraint(trait, labrid$tree,
                                          regime=two_shifts, 
                                          alpha=(ou3@sqrt.alpha)^2,
                                          sigma=ou3@sigma) 




  loglik(wright_twoshifts)-loglik(wright_intra)


  results <- matrix(NA, nrow=4, ncol=3, dimnames = list(c("ouch", "brownie", "release", "wright"), c("phar", "intra", "twoshifts")))

  results[1,1] <- loglik(ouch_phar)
	results[2,1] <- loglik(brownie_phar)
  results[3,1] <- loglik(release_phar)
	results[4,1] <- loglik(wright_phar)

  results[1,2] <- loglik(ouch_intra )
	results[2,2] <- loglik(brownie_intra )
  results[3,2] <- loglik(release_intra)
	results[4,2] <- loglik(wright_intra )

  results[1,3] <- loglik(ouch_twoshifts )
	results[2,3] <- loglik(brownie_twoshifts)
  results[3,3] <- loglik(release_twoshifts)
	results[4,3] <- loglik(wright_twoshifts)


#  barplot(results, xlim=c(0, 80), col=c("thistle", "khaki", "pink","palegreen"), horiz=TRUE, beside=TRUE)
  social_plot(barplot(t(results), xlim=c(0, 80), col=c("thistle", "khaki", "palegreen"), horiz=TRUE, beside=TRUE, main=trait_name), tag=tag, comment=trait_name)



# specify a regime
# function of a model m
  n_regimes <- length(levels(m$regimes))
  n_traits <- length(getParameters(m)) # -1 # -1 if one is for convergence
  parameters <- matrix(NA, ncol=n_traits, nrow=n_regimes)
  parameters[i,j] <- getParameters(brownie_phar) 


}

cpu <- 9
nboot <- 160
sfInit(parallel=TRUE, cpu=cpu)
sfExportAll()
sfLibrary(wrightscape)  # need all this just to export wrightscape?
sfLibrary(pmc)
sfLibrary(socialR)


sfSapply(1:9, function(i) try(dummy(i)) )



 # Can we do better?
# wright_twoshifts_ <- wright(trait, labrid$tree, regime=two_shifts, alpha=c(wright_intra$alpha, wright_intra$alpha[2]), sigma=c(wright_intra$sigma, wright_intra$sigma[2]))


## Doing the likelihood optimization in C instead.  needs robustness testing, convergence conditions still
#  ws2_phar <- wrightscape(trait, labrid$tree, regime=labrid$regimes, (ou2_phar@sqrt.alpha)^2, ou2_phar@sigma, theta=ou2_phar@theta[[1]])
#  ws2_intra <- wrightscape(trait, labrid$tree, regime=intramandibular, (ou2_intra@sqrt.alpha)^2, ou2_intra@sigma, theta=ou2_intra@theta[[1]])
#  ws2_twoshifts <- wrightscape(trait, labrid$tree, regime=intramandibular, (ou3@sqrt.alpha)^2, ou3@sigma, theta=ou3@theta[[1]])


## Make this pretty using pmc and ape plot tools
plt <- function(){
  phylo_phar <- convert(ou2_phar)
  phylo_intra <- convert(ou2_intra)
  phylo_twoshifts <- convert(ou3)
  ou3 <- hansen(trait, labrid$tree, regime=two_shifts, ou2_phar@sqrt.alpha, sigma=ou2_phar@sigma )
  par(mfrow=c(1,3))
  plot(phylo_phar, edge.color = treepalette(phylo_phar), edge.width=2, cex=1.0, show.tip.label=FALSE)
  plot(phylo_intra, edge.color = treepalette(phylo_intra), edge.width=2, cex=1.0, show.tip.label=FALSE)
  plot(phylo_twoshifts, edge.color = treepalette(phylo_twoshifts), edge.width=2, cex=1.0, show.tip.label=FALSE)
}
#social_plot(plt(), file="labrids.png", tag=tag)


#plot(labrid$tree, regimes=two_shifts)


## Do some modelchoice
#out <- montecarlotest(brownie_phar, release_phar, cpu=cpu,nboot=nboot, GetPar=F) 
#social_plot(plot(out), tag="phylogenetics wrightscape labrids", comment="brownie vs release on pharyngeal shift pt, trait = gape")
#out2 <- montecarlotest(release_phar, release_intra, cpu=cpu,nboot=nboot, GetPar=F) 
#social_plot(plot(out2), tag="phylogenetics wrightscape labrids", comment="release on pharyngeal vs intramandibular shift, trait=gape")
#out3 <- montecarlotest(release_intra, release_twoshifts, cpu=cpu,nboot=nboot, GetPar=F) 
#social_plot(plot(out3), tag="phylogenetics wrightscape labrids", comment="release on pharyngeal vs intramandibular shift, trait=gape")



#})



