# loop_models_traits_regimes.R

fit <- function(model, traits, regimes, tree, alpha=0.1, sigma=0.1){
  switch(model,
   hansen =   hansen(traits, tree, regimes, alpha, sigma),
   ouch =       ouch(traits, tree, regimes, alpha, sigma),
   brownie = brownie(traits, tree, regimes, sigma=sigma),
   wright =   wright(traits, tree, regimes, alpha=alpha, sigma=sigma),
   release_constraint = 
  release_constraint(traits, tree, regimes, alpha=alpha, sigma=sigma))
}



results <- function(models, traits, regimes, tree){ 
  lapply(1:length(traits), function(j){
    lapply(1:length(regimes), function(k){
      ## The loop over models will use a for loop so it 
      ## can try alpha/sgima values from earlier models
      out <- vector("list", length=length(models))
      for(i in 1:length(models)){
        
        ##  Hackish treatment of models that don't have regimes.  
        if (models[[i]]=="brown"){
            out[[i]] <- brown(traits[j], tree)
        } else if (models[[i]]=="ou1"){
          out[[i]] <- hansen(traits[j], tree,  regime=labrid$noregimes,
                             .01, .01)
        } else if(models[[i]]=="hansen"){
          out[[i]] <- try( fit(models[[i]], traits[j], regimes[[k]],
                               tree=tree) )
        } else {
          print(paste("model = ", models[[i]], "trait = ", names(traits[j]),
                regime=names(regimes[[k]])))
          out[[i]] <- try( fit(models[[i]], traits[j], regimes[[k]],
                               tree=tree, alpha=(out[[2]]@sqrt.alpha)^2,
                               sigma=out[[2]]@sigma))

          ## attempt previous result as starting conditions
          if(is(out[[i]], "try-error")){
          }
        }
      }
      out
    })
  })
}


