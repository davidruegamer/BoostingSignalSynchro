#######################################################################################
### Simulation code for comparison of different parametrizations / model w/wo shrinkage
#######################################################################################

source("0_libs_funs.R")

######### settings
n <- c(80, 320)
obsPerTra <- c(60)
nrFacLev <- c(4,10)
SNR <- c(1/10, 1, 10)
setup = c("histGameIA")
nrSims = 100
nrBS = 25

######### generate all combinations of different settings
setupDF <- expand.grid(list(nrFacLev=nrFacLev,
                            obsPerTra=obsPerTra,
                            n=n,
                            SNR=SNR))

resSim <- vector("list",nrow(setupDF))

# do for all settings
for(i in length(resSim):10){
  
  ######### extract settings
  nrFacLev = setupDF$nrFacLev[i]
  obsPerTra = setupDF$obsPerTra[i]
  n = setupDF$n[i]
  SNR = setupDF$SNR[i]
  
  ######### generate data
  dat <- dataGenProc(n = n,
                     obsPerTra = obsPerTra,
                     SNR = SNR,
                     seed = 111,
                     setup = setup,
                     nrOfResp = nrSims,
                     nrFacEf = nrFacLev,
                     nrRanEf = 8
  )
  
  foldMat <- apply(sapply(1:nrBS, function(y) c(sapply(1:nrFacLev, 
                             function(x) 
                               sample(1:(n/nrFacLev), replace = T))) + 
                      rep(seq(0,n,by=(n/nrFacLev))[-(nrFacLev+1)], each=(n/nrFacLev))),2,
                   function(x)rep(c(table(factor(x, levels=1:n))),obsPerTra))
  
  ######### parallelize over repetitions
  simDF <- mclapply(1:nrSims,function(nrSim){
    
    tryCatch({
    ######### use the ith response matrix
    dat$Yi <- dat$Y[[nrSim]]
    dat$Yvec <- Yvec <- as.vector(dat$Yi)
    
    # fit models
    ###############################################################
    
    ######### fit model 1 (with deviations from a main historical effect)
    modCurrent <- FDboost(Yi ~ 1 + bhistx(X1h, df=45, knots=5, differences=2, standard='length') + 
                            bhistx(X1h, df=45, knots=5, differences=2, standard='length') %X% 
                            bolsc(g2, index=repIDx, df=1) + 
                            bolsc(g2, df=nrFacLev),
                          timeformula = ~ bbs(t, df=(45/nrFacLev)), data=dat,
                          control=boost_control(mstop=700,nu=0.1)
    )
    
    ######### validate
    cvr1 <- cvrisk(modCurrent, folds = foldMat, mc.cores=25)
    modCurrentFin <- modCurrent[mstop(cvr1)]
    
    ######### fit model 2 (with deviations from a main historical effect but with no shrinkage for g2)
    modWoShrink <- FDboost(Yi ~ 1 + bhistx(X1h, df=45, knots=5, differences=2, standard='length') + 
                             bhistx(X1h, df=45, knots=5, differences=2, standard='length') %X% 
                             bolsc(g2, index=repIDx, df=1, 
                                   K = matrix(0, ncol = length(levels(dat$g2)), nrow = length(levels(dat$g2)))) + 
                             bolsc(g2, df=nrFacLev) %A% bbs(t, df=(45/nrFacLev)),
                           timeformula = ~ bbs(t, df=(45/nrFacLev)), data=dat,
                           control=boost_control(mstop=700,nu=0.1)
    )
    
    ######### validate
    cvr2 <- cvrisk(modWoShrink, folds = foldMat, mc.cores=25)
    modWoShrinkFin <- modWoShrink[mstop(cvr2)]
    
    ######### fit model 3 (with deviations from a main historical effect but no shrinkage trough pen.mat.)
    modXa0 <- FDboost(Yi ~ 1 + bhistx(X1h, df=45, knots=5, differences=2, standard='length') + 
                        bhistx(X1h, df=(45/(nrFacLev-1)), 
                               knots=5, differences=2, standard='length') %Xa0% 
                        bolsc(g2, index=repIDx, df=nrFacLev-1) + 
                        bolsc(g2, df=nrFacLev) %A% bbs(t, df=(45/nrFacLev)),
                      timeformula = ~ bbs(t, df=45/nrFacLev), data=dat,
                      control=boost_control(mstop=700,nu=0.1)
    )
    
    ######### validate
    cvr3 <- cvrisk(modXa0, folds = foldMat, mc.cores=25)
    modXa0Fin <- modXa0[mstop(cvr3)]

    
    ######### fit model 4 (only with game-condition specific effect)
    modWoC <- FDboost(Yi ~ 1 + bhistx(X1h, df=45, knots=5, differences=2, standard='length') %X% 
                        bols(g2, index=repIDx, df=1, contrasts.arg = "contr.dummy") + 
                        bolsc(g2, df=nrFacLev),
                      timeformula = ~ bbs(t, df=45/nrFacLev), data=dat,
                      control=boost_control(mstop=700,nu=0.1)
    )
    
    ######### validate
    cvr4 <- cvrisk(modWoC, folds = foldMat, mc.cores=25)
    modWoCFin <- modWoC[mstop(cvr4)]

    
    ######### fit model 5 (only with game-condition specific effect and no penalty)
    modXa0woC <- FDboost(Yi ~ 1 + bhistx(X1h, df=45/nrFacLev, knots=5, 
                                         differences=2, standard='length') %Xa0% 
                           bols(g2, index=repIDx, df=nrFacLev, contrasts.arg = "contr.dummy") + 
                           bolsc(g2, df=nrFacLev) %A% bbs(t, df=(45/nrFacLev)),
                         timeformula = ~ bbs(t, df=45/nrFacLev), data=dat,
                         control=boost_control(mstop=700,nu=0.1)
    )
    
    ######### validate
    cvr5 <- cvrisk(modXa0woC, folds = foldMat, mc.cores=25)
    modXa0woCFin <- modXa0woC[mstop(cvr5)]
    
    ###########################################################################
    
    mods <- list(modCurrentFin,
                 modWoShrinkFin,
                 modXa0Fin,
                 modWoCFin,
                 modXa0woCFin
                 )
    
    reliMsesPerMod <- unlist(mclapply(1:length(mods), function(j){
      
      modFin <- mods[[j]]
      
      ######### get relimses
      
      relimse = vector("list",nrFacLev)
      
      
      ## Main Effect
      
      relimse <- mean(sapply(1:nrFacLev,function(fl){
        
        trueX1eff <- dat$trueEffHist(dat$s,dat$t) + 
          dat$trueEffHistGame(dat$s,dat$t)*dat$trueEffHistGameFac[fl]
        
        addOne <- !j%in%4:6
        
        predEff <- try(coef(modFin, which=2:(2+as.numeric(addOne)), n1=obsPerTra, n2=obsPerTra)$smterms)
        if(class(predEff)=="try-error") return(NA)
        
        predEff1 <- if(addOne) predEff[[1]]$value + predEff[[2]][[fl]]$value else predEff[[1]][[fl]]$value
        
        trueX1eff[lower.tri(trueX1eff)] <- NA
        predEff1[lower.tri(predEff1)] <- NA
        
        sum(c(((predEff1-trueX1eff)^2)),na.rm = T)/sum(c(trueX1eff^2),na.rm = T)
        
      }))
      
      return(relimse) # averaged over factor levels
      
    }, mc.cores=5))
    
    return(cbind(data.frame(relimse = reliMsesPerMod,
                            mstopIter = sapply(mods,mstop), 
                            modType = c("modCurrentFin",
                                        "modWoShrinkFin",
                                        "modXa0Fin",
                                        "modWoCFin",
                                        "modXa0woCFin"),
                            nrSim = nrSim), setupDF[i,]))
    
    }, error=function(e)return(e))
    
  }, mc.cores=1)
  
  saveRDS(simDF, file=paste0("results/perfParam/","simNr",i,".RDS"))
  resSim[[i]] <- simDF
  
}

######### save results
saveRDS(resSim,file="iaG_comp_param.RDS")