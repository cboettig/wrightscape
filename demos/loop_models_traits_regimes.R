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



fit_all <- function(models, traits, regimes, tree){ 
  fits <- lapply(1:length(traits), function(j){
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


        ## Simple hansen model
        } else if(models[[i]]=="hansen"){
          out[[i]] <- try( fit(models[[i]], traits[j], regimes[[k]],
                               tree=tree) )

        ##   
        } else {
          ## Reporting, optional
#          print(paste("model = ", models[[i]], "trait = ", names(traits[j]),
#                      "regime=", names(regimes[k])))
          ## use hansen to start with good parameters 
          hansen <- try(fit("hansen", traits[j], regimes[[k]], tree=tree))
          ## hansen will sometimes give very large or negative 

          ## Fit one of the generalized models using the initial guess from hansen
          out[[i]] <- try( fit(models[[i]], traits[j], regimes[[k]],
                               tree=tree, alpha=(min(10, hansen@sqrt.alpha^2)),
                               sigma=hansen@sigma))
          ## If errors, attempt default starting conditions 
          if(is(out[[i]], "try-error")){
            out[[i]] <- try( fit(models[[i]], traits[j], regimes[[k]],
                               tree=tree, alpha=0.01,
                               sigma=0.01))
          }
        }
      }
      out  
    })
  })
  list(fits=fits, models=models, traits=traits, regimes=regimes, tree=tree)
}


llik_matrix <- function(results){

  models <- results$models
  regimes <- names(results$regimes)
  traits <- names(results$traits)

  all <- vector("list", length(traits))
  for(j in 1:length(traits)){
    M <- matrix(NA, nrow=length(regimes), ncol=length(models)) 
    for(i in 1:length(models)){
      for(k in 1:length(regimes)){
        a <- results$fit[[j]][[k]][[i]]
        if (is(a, "try-error")){
          M[k,i] <- NA
        } else if(is(a, "hansentree") | is(a, "browntree")){
          M[k,i] <- a@loglik
        } else if(is(a, "multiOU")){
          M[k,i] <- loglik(a)
        }
      }
      rownames(M) <- regimes
      colnames(M) <- models
    }
    all[[j]] <- M
  }
  names(all) <- traits
  all
}


alpha_matrix <- function(results, trait_id, k = 1){
  k <- 1
  j <- trait_id 
  models <- results$models
  paintings <- names(results$regimes)
  traits <- names(results$traits)

  groups <- levels(results$regimes[[k]])
  n_groups <- length(groups)

  M <- matrix(NA, nrow=n_groups, ncol=length(models)) 
  for(i in 1:length(models)){
      a <- results$fit[[j]][[k]][[i]]
    for(m in 1:n_groups){
      if (is(a, "try-error") | is(a, "browntree")){
        print(class(a))
        M[m,i] <- NA
      } else if(is(a, "hansentree")){
        M[m,i] <- (a@sqrt.alpha)^2
      } else if(is(a, "multiOU")) {
        if(a$submodel != "brownie"){
        M[m,i] <- a$alpha[m]
        } else {
         M[m,i] <- NA
        }
      }
    }
    colnames(M) <- models
    rownames(M) <- groups
  }
  M
}


