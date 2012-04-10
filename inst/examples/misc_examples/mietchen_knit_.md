# mietchen.R

 load some useful libraries and the phylogeny
``` {r }
require(auteur)
data(primates)
````


load in Daniel's data
``` {r }
dat <- read.csv("../../../data/Primate_brain_comparisons.csv")
````
Concatenate Genus and species names
``` {r }
SpeciesNames <- sapply(1:length(dat[[1]]), 
  function(i) 
  paste(dat[[1]][i], "_", dat[[2]][i], sep=""))
````
 seems these two genus names differ in spelling in the data and tree,
 let's just fix them manually


``` {r }
SpeciesNames <- gsub("Saguinas", "Saguinus", SpeciesNames)
SpeciesNames <- gsub("Presbytus", "Presbytis", SpeciesNames)
````

Name the data rows as "Genus_species", matching the tree tip label convention


``` {r }
rownames(dat) <- SpeciesNames
````


get just the quantitative trait data

``` {r }
dat <- dat[3:7]
````

 drop all tips that don't have data, or data that doesn't have tips in tree

``` {r }
primate_data <- treedata(primates$phy, dat)
dat <- primate_data$data
phy <- primate_data$phy
```` 
Classic independent contrasts for a phylogenetic correction to the 
 estimate of correlations in brain weight with body size

``` {r }
x <- pic(log(dat[,"Body_weight"]), phy)
y <- pic(log(dat[,"Brain_weight"]), phy)
summary(lm(y~x-1))
````


 get all ancestral states. Note that Ancestral state estimates
 are highly uncertain and generally distrusted, see Schluter et al 1997.
``` {r }
ancestral_states <- lapply(1:dim(dat)[2], function(i) ace(dat[,i], phy))
````


 get the BM diversification rates for each trait
``` {r }
diversification_rates <- fitContinuous(phy, dat)
````


## Estimate a shift in diversification rate using AUTEUR           
run two short reversible-jump Markov chains
(create some random strings for temporary file names)
``` {r }
r=paste(sample(letters,9,replace=TRUE),collapse="")
````

 run four short MCMC chains to search for a change point in brain weight
``` {r }
require(snowfall)
sfInit(parallel=TRUE, cpu=4)
sfLibrary(auteur)
sfExportAll()
out <- sfLapply(1:4, 
         function(x) rjmcmc.bm(phy=phy, dat=dat[,"log_brain.weight"],
          ngen=100000, sample.freq=10, prob.mergesplit=0.1, simplestart=TRUE,
          prop.width=1, fileBase=paste(r,x,sep=".")))
````


collect directories

``` {r }
dirs=dir("./",pattern=paste("BM",r,sep="."))
pool.rjmcmcsamples(base.dirs=dirs, lab=r)
````
view contents of .rda
``` {r }
 load(paste(paste(r,"combined.rjmcmc",sep="."), 
      paste(r,"posteriorsamples.rda",sep="."),sep="/"))
 print(head(posteriorsamples$rates))
 print(head(posteriorsamples$rate.shifts))
````


plot Markov sampled rates

``` {r }
  shifts.plot(phy=phy, base.dir=paste(r,"combined.rjmcmc",sep="."), burnin=0.5, legend=TRUE, edge.width=4, x.lim = c(0,60))
````





