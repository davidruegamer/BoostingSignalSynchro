### function to perform bootstrap with appropriate resampling
### 
### if retBootCoefs = FALSE, the function will process the
### bootstrap results in order to obtain the quantiles
### for coefficients specified by which
bootMod <- function(mod, 
                    nrBoots, # number of bootstrap iterations
                    vars, argvals, idvars, # variables passed to reweightData (see ?reweightData)
                    Yname = "Yi", # name of the response
                    by = NULL, # for sampling on specific levels
                    data,
                    quantiles = list(c(0.005,0.995),
                                     c(0.025,0.975),
                                     c(0.05,0.95),
                                     c(0.1,0.9)), 
                    which = 2, # which coefficient should be extracted
                    obsPerTra = 40, # number of observations per trajectory for plot function
                    indFun = function(x)x$value, # function which extracts the coefficients
                    mcCores = 20, # cores to use
                    retBootCoefs = TRUE, 
                    ... # further arguments passed to cvrisk
                    )
{
  
  obsPerTra <- min(obsPerTra, 80) # due to restrictions in coef.FDboost
  
  ###############################################################################
  ###############################################################################
  # do the following for all nrBoots iterations
  bootCoefs <- mclapply(1:nrBoots, function(i){
    
    set.seed(i)
    
    # if by is given, sampling is done on data[[by]] levels
    if(!is.null(by)){
      
      newLev <- sample(levels(data[[by]]), replace=T)
      ind <- unlist(lapply(newLev, function(nl) which(data[[by]] == nl)))
      
    }else{ # ordinary sampling
      
      ind <- sample(1:mod$ydim[1], replace = T)
      
    }
    
    # perform resampling using the index ind
    datNew <- reweightData(data, vars = vars, argvals = argvals, 
                           index = ind, idvars = idvars)
    
    # if meaningful, drop empty levels of by variable
    if(!is.null(by)) datNew[[by]] <- droplevels(datNew[[by]])
    
    # refit the model using the new data
    modN <- update(mod, data = datNew)
    # extract modstop
    mstopModN <- mstop(modN)
    
    # create folds
    ppmat <- if(!is.null(by)) createRandomRespFolds(ranVar = datNew[[by]], sLength = modN$ydim[2]) else
      cvLong(modN$id, weights = rep(1,l=length(modN$id)), type = "kfold")
    
    # calculat mstop via cross-validation
    cvr <- cvrisk(modN, folds = ppmat, grid = 1:mstopModN, papply = lapply, 
                  fun = NULL, corrected = TRUE, ...)
    
    # extract mstop
    mcvr <- mstop(cvr)
    
    # define final model
    modN <- modN[mcvr]
    
    # extract coefficient path of final model
    selCourse <- selected(modN)
    # check if which was selected 
    which <- intersect(selCourse, which)
    # get coefficients
    coefsList <- lapply(which, function(w) coef(modN, which = w, n1 = obsPerTra, n2 = obsPerTra)$smterms[[1]])
    # return coefficients and which
    return(list(cl=coefsList,w=which))
    
  },mc.cores=mcCores)
  ###############################################################################
  ###############################################################################

  # return coefficients 
  if(retBootCoefs) return(bootCoefs)
  
  # else calculate quantiles of coefficients
  # -> recommendation: do this manually as this is not computationally costly
  # and it's much easier to intervene when things go wrong
  
  #############################################################################################
  #############################################################################################
  #############################################################################################
  
  wt <- lapply(bootCoefs, "[[", "w")
  which <- Reduce("unique", wt[sapply(wt, length)>0])
  bootCoefs <- lapply(bootCoefs, "[[", "cl")
  
  retList <- vector("list", length(which))
  i <- 1
  namesEffects <- sapply(bootCoefs[[which(sapply(bootCoefs, length) == length(which))[1]]], 
                         "[[", "main")
  
  for(w in which){
    
    nameEffequalsWhich <- namesEffects[i]
    
    flatCoefs <- lapply(bootCoefs, function(bc){ 
      
      effInBc <- sapply(bc,"[[","main")
      
      if(!nameEffequalsWhich%in%effInBc) return(NULL)
      
      indEff <- which(effInBc==nameEffequalsWhich)

      if(!is.null(bc[[indEff]]$numberLevels)){
        return(lapply(1:bc[[indEff]]$numberLevels,function(j)c(indFun(bc[[indEff]][[j]]))))
      }else{
        return(c(indFun(bc[[indEff]])))
      }
      
    })
    
    if(is.list(flatCoefs)) flatCoefs <- flatCoefs[!sapply(flatCoefs,is.null)]
    if(length(flatCoefs)==0){ 
      warning(paste0("Effect ",w," was never selected")) 
      retList[[i]] <- NULL 
      }
    
    indLev <- !is.null(nrLev <- bootCoefs[[which.max(sapply(bootCoefs,length))]][[i]]$numberLevels)
    
    ncolD <- mod$ydim[2]
    if(!indLev & is.list(flatCoefs)) flatCoefs <- do.call("cbind",flatCoefs)
    
    if(indLev){
      
        if(is.list(quantiles)){
          
          
          bootCIs <- nonZero <- vector("list",length(quantiles))
          
          for(ql in 1:length(quantiles)){
            
            qCoefs <- lapply(1:nrLev, function(lev)apply(do.call("rbind",lapply(flatCoefs,"[[",lev)),2,quantile,
                                                         quantiles[[ql]]))
            bootCIs[[ql]] <- lapply(1:nrLev, function(lev) list(q1=matrix(qCoefs[[lev]][1,],ncol=obsPerTra),
                                                          q2=matrix(qCoefs[[lev]][2,],ncol=obsPerTra)))
            nonZero[[ql]] <- lapply(1:nrLev, function(lev) Reduce("*",bootCIs[[ql]][[lev]])>0)
            
          }
          
        }else{
          
          qCoefs <- lapply(1:nrLev, function(lev)apply(do.call("rbind",lapply(flatCoefs,"[[",lev)),2,
                                                       quantile,quantiles))
          bootCIs <- lapply(1:nrLev, function(lev) list(q1=matrix(qCoefs[[lev]][1,],ncol=obsPerTra),
                                                        q2=matrix(qCoefs[[lev]][2,],ncol=obsPerTra)))
          nonZero <- lapply(1:nrLev, function(lev) Reduce("*",bootCIs[[lev]])>0)
          
        }

      
    }else{

        if(dim(flatCoefs)[1]>dim(flatCoefs)[2]) flatCoefs <- t(flatCoefs)
        
        if(is.list(quantiles)){
          
          bootCIs <- nonZero <- vector("list",length(quantiles))
          
          for(ql in 1:length(quantiles)){
            
            qCoefs <- apply(flatCoefs,2,quantile,quantiles[[ql]])
            bootCIs[[ql]] <- list(q1=matrix(qCoefs[1,],ncol=obsPerTra),
                            q2=matrix(qCoefs[2,],ncol=obsPerTra))
            nonZero[[ql]] <- Reduce("*",bootCIs[[ql]])>0
            
          }
          
        }else{
          
          qCoefs <- apply(flatCoefs,2,quantile,quantiles)
          bootCIs <- list(q1=matrix(qCoefs[1,],ncol=obsPerTra),
                          q2=matrix(qCoefs[2,],ncol=obsPerTra))
          nonZero <- Reduce("*",bootCIs)>0
        
        }

    }

    retList[[i]] <- list(bootCIs = bootCIs, 
                         nonZero = nonZero)
    i <- i+1
    
  }

  retList <- retList[!sapply(retList,is.null)]
  
  return(retList)
  
}