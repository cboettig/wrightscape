# loop_models_traits_regimes.R

fit <- function(model, traits, regimes, tree, alpha=0.1, sigma=0.1, ...){
  switch(model,
   hansen =   hansen(traits, tree, regimes, alpha, sigma, ...),
   ouch =       ouch(traits, tree, regimes, alpha, sigma, ...),
   brownie = brownie(traits, tree, regimes, sigma=sigma, ...),
   wright =   wright(traits, tree, regimes, alpha=alpha, sigma=sigma, ...),
   release_constraint = 
  release_constraint(traits, tree, regimes, alpha=alpha, sigma=sigma, ...))
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
                             .01, .01, maxit=5000)


        ## Simple hansen model
        } else if(models[[i]]=="hansen"){
          out[[i]] <- try( fit(models[[i]], traits[j], regimes[[k]],
                               tree=tree, .01, .01, control=list(maxit=5000)) )

        ##   
        } else {
          ## Reporting, optional
#          print(paste("model = ", models[[i]], "trait = ", names(traits[j]),
#                      "regime=", names(regimes[k])))


          ## use hansen to start with good parameters 
          hansen <- try(fit("hansen", traits[j], regimes[[k]], tree=tree,
                            maxit=5000))


          ## Fit one of the generalized models using the initial guess from hansen
          out[[i]] <- try( fit(models[[i]], traits[j], regimes[[k]],
                               tree=tree, alpha=c(5,.01),
#(min(10, hansen@sqrt.alpha^2)),
                               sigma=hansen@sigma, method ="SANN",
                               control=list(maxit=50000,temp=15,tmax=20)))

          ## If errors, attempt default starting conditions 
          if(is(out[[i]], "try-error")){
            warning("first attempt to fit failed, trying another routine")
            out[[i]] <- try( fit(models[[i]], traits[j], regimes[[k]],
                               tree=tree, alpha=0.01,
                               sigma=0.01, control=list(maxit=5000)))
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



sort_lik <- function(likmat, model_names = c("bm", "ou", "theta",
                      "sigma", "gen", "alpha")){
## Takes output of llik_matrix and:
## Reorganize the data to list by ascending score
## Rename models (with short, parameter-based names) 
## rescale the data relative to weakest model 
lliks <- vector("list", length=length(likmat))
for(i in 1:length(likmat)){
  tmp <- likmat[[i]]
  names(tmp) <- model_names
  tmp <- sort(tmp)
  shift <- tmp[1] 
  for(j in 1:length(tmp)){
    print(tmp[j])
    tmp[j] <- tmp[j]-shift
    print(tmp[j])
  }
  lliks[[i]] <- tmp
  print(tmp)
}
names(lliks) <- names(likmat)
lliks
}


alpha_traits <- function(results, model_name="release_constraint"){
# Takes the output of fit_all and grabs alpha values for specified model
  k <- 1 # the regime index
  groups <- levels(results$regimes[[k]])
  n_groups <- length(groups)
  i <- match(model_name, results$models)
  M <- matrix(NA, ncol=n_groups, nrow=length(results$traits))
  for(j in 1:length(results$traits)){
    a <- results$fit[[j]][[k]][[i]]
    M[j,] <- a$alpha
  }
  rownames(M) <- names(results$traits)
  colnames(M) <- groups
  M
}


conv<- function(results){
## Takes the output of fit_all and reports which ones didn't converge
  for(j in 1:length(results$traits)){
    for(k in 1:length(results$regimes)){
      for(i in 1:length(results$models)){
        a <- results$fit[[j]][[k]][[i]]
        if (is(a, "multiOU"))
          if (a$convergence != 0)
            print(paste("trait ", names(results$traits)[j], 
            ", regime ", names(results$regimes)[k], ", model ", 
            results$models[i], " didn't converge", sep=""))
      }
    }
  }
}





alpha_matrix <- function(results, trait_id, k = 1){
## Takes the output of fit_all and computs the matrix of alphas 
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



