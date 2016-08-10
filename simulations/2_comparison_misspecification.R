############################################################################
### Simulation code for section 4.3 (missspecification)
############################################################################

source("0_libs_funs.R", chdir = T)
if(length(list.files("results")) == 0) dir.create("results")

nrSims = 100
## if you just want to test the code:
if(FALSE) nrSims = 2

### core usage
coresCV = 1
coresSettings = 9


######### settings
addME <- FALSE # TRUE
nuC <- c(0.1)
n <- c(80, 160, 320)
obsPerTra <- c(#20, 
  40#, 60
)
SNR <- c(0.1, 1, 10)
nrRanEf <- c(10)
setupMiss = data.frame(dgp = c(#"histGameIA", "histGameIA", 
                               #"histRandIA", "histRandIA", 
                               "full"),
                       fit = c(#"histOnly", "histAndGame", 
                               #"histOnly", "histAndRand", 
                               "withoutDoubleVar"))
setupComb <- 1:nrow(setupMiss)

######### generate all combinations of different settings
setupDF <- expand.grid(list(setupNr=setupComb,
                            n=n,
                            SNR=SNR#,
                            #obsPerTra=obsPerTra,
                            #nrRanEf=nrRanEf
))
# setupDF$setup <- as.character(setupDF$setup)

######### parallelize over different settings
resSim <- mclapply(1:nrow(setupDF),function(i){
  
  ######### extract settings
  setupDGP = as.character(setupMiss[setupDF$setupNr[i],1])
  setupFIT = as.character(setupMiss[setupDF$setupNr[i],2])
  n = setupDF$n[i]
  SNR = setupDF$SNR[i]
  #obsPerTra = setupDF$obsPerTra[i]
  #nrRanEf = setupDF$nrRanEf[i]
  
  ######### generate data
  dat <- dataGenProc(n = n,
                     obsPerTra = obsPerTra,
                     SNR = SNR,
                     seed = 12,
                     setup = setupDGP,
                     nrOfResp = nrSims,
                     nrRanEf = nrRanEf,
                     nrFacEf = 4
  )
  
  ######### get model specification
  baseForm <- "Yi ~ 1" 
  terms <- c("bolsc(g2, df=6)",   #1
             "brandomc(g3, df=6)", #2
             "bols(g2, df=2) %Xc% brandom(g3, df=3)",  #3
             "bhistx(X1h, df=15, knots=5, differences=2, standard='length')", #4
             "bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% bolsc(g2, index=repIDx, df=1)",
             "bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% brandomc(g3, index=repIDx, df=1)",
             "bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% myBlg"
  )
  
  ind <- switch(setupFIT,
                histOnly = 4,
                histAndGame = c(4,1),
                histAndRand = c(4,2),
                histGameRand = c(4,1,2,3),
                histGameIA = c(4,5,1),
                histRandIA = c(4,6,2),
                full = c(4,5,6,7,1,2,3),
                withoutDoubleVar = c(4,5,6,1,2,3)
  )
  
  indTru <- switch(setupDGP,
                histOnly = 4,
                histAndGame = c(4,1),
                histAndRand = c(4,2),
                histGameRand = c(4,1,2,3),
                histGameIA = c(4,5,1),
                histRandIA = c(4,6,2),
                full = c(4,5,6,7,1,2,3),
                withoutDoubleVar = c(4,5,6,1,2,3)
  )
  
  fff <- as.formula(paste(baseForm,paste(terms[ind],collapse=" + "),sep=" + "))
  fffTru <- as.formula(paste(baseForm,paste(terms[indTru],collapse=" + "),sep=" + "))
  
  simDF <- vector("list",nrSims)
  
  ######### do for nrSims repetitions
  for(nrSim in 1:nrSims){
    
    ######### use the ith response matrix
    dat$Yi <- if(nrSims>1) dat$Y[[nrSim]] else dat$Y
    dat$Yvec <- Yvec <- as.vector(dat$Yi)
    
    myBlg <<-  with(dat,(bols(g2, index = repIDx, df = 1) %Xc%
                           brandom(g3, index = repIDx, df = 1)))
    
    
    
    ########### Fit misspecified and true model
    
    mod1 <- FDboost(fff,
                    timeformula = ~ bbs(t, df = 2.5), 
                    data = dat,
                    control = boost_control(mstop = 2500, nu = 0.1)
    )
    
    gridEnd = 2500
    gridStart = 1
    
    ######### generate appropriate folds
    ppmat <- createRandomRespFolds(ranVar = dat$g3, folds = 10, sLength = obsPerTra)
    
    ######### validate
    cvr <- if(setupFIT %in% c("histAndRand", "withoutDoubleVar")) 
      cvrisk(mod1, grid = gridStart:gridEnd, 
             folds = ppmat, mc.cores = coresCV) else
               cvrisk(mod1, folds = cvLong(id = mod1$id, 
                                           weights = model.weights(mod1),
                                           B = 10))
    
    modFin1 <- mod1[mstop(cvr)]
    
    
    mod2 <- FDboost(fffTru,
                    timeformula = ~ bbs(t, df = 2.5), 
                    data = dat,
                    control = boost_control(mstop = 2500, nu = 0.1)
    )
    
    gridEnd = 2500
    gridStart = 1
    
    ######### generate appropriate folds
    ppmat <- createRandomRespFolds(ranVar = dat$g3, folds = 10, sLength = obsPerTra)
    
    ######### validate
    cvr <- if(setupDGP %in% c("histAndRand", "withoutDoubleVar")) 
      cvrisk(mod2, grid = gridStart:gridEnd, 
             folds = ppmat, mc.cores = coresCV) else
               cvrisk(mod2, folds = cvLong(id = mod2$id, 
                                           weights = model.weights(mod2),
                                           B = 10))
                    
    modFin2 <- mod2[mstop(cvr)]
    
    
    ######### compare models
    
    findEffects <- which(c(4:7)%in%ind)
    selCourse <- selected(modFin1)
    selCourseTru <- selected(modFin2)
    
    relimseMain <- NA
    relimseIAGame <- as.list(rep(NA,length(levels(dat$g2))))
    relimseIARan <- as.list(rep(NA,length(levels(dat$g3)))) 
    
    relimseMainTru <- NA
    relimseIAGameTru <- as.list(rep(NA,length(levels(dat$g2))))
    relimseIARanTru <- as.list(rep(NA,length(levels(dat$g3)))) 
    
    ######### get relimses
    
    # main effect is fitted in every setting
    
    if(1%in%findEffects & 2%in%selCourse){
      
      ## Main Effect
      
      trueX1eff <- dat$trueEffHist(dat$s,dat$t)
      trueX1eff[trueX1eff==0] <- NA
      predEff <- coef(modFin1,which=2,n1=obsPerTra,n2=obsPerTra)$smterms[[1]]$value
      predEff[predEff==0] <- NA
      
      relimseMain <- sum(c(((predEff-trueX1eff)^2)),na.rm=T)/sum(c(trueX1eff^2),na.rm = T)
      
      
    }
    
    if(1%in%findEffects & 2%in%selCourseTru){
      
      ## Main Effect
      
      trueX1eff <- dat$trueEffHist(dat$s,dat$t)
      trueX1eff[trueX1eff==0] <- NA
      predEff <- coef(modFin2,which=2,n1=obsPerTra,n2=obsPerTra)$smterms[[1]]$value
      predEff[predEff==0] <- NA
      
      relimseMainTru <- sum(c(((predEff-trueX1eff)^2)),na.rm=T)/sum(c(trueX1eff^2),na.rm = T)
      
      
    }
    
    if(setupFIT == "withoutDoubleVar"){
    
      # Interaction Effect with Categorical Covariate
      
      if(2%in%findEffects & 3%in%selCourse){
        
        ccc <- coef(modFin1,which=3,n1=obsPerTra,n2=obsPerTra)$smterms
        
        for(nr in 1:length(levels(dat$g2))){
          
          truth1 <- dat$trueEffHistGame(dat$s,dat$t)*dat$trueEffHistGameFac[nr]
          predEff1 <- ccc[[1]][[nr]]$value
          predEff1[predEff1==0] <- NA
          truth1[lower.tri(truth1)] <- NA
        
          relimseIAGame[[nr]] <- sum(c(((predEff1-truth1)^2)),na.rm=T) / sum(c(truth1^2),na.rm = T)
          
        }
        
      }
      
      if(2%in%findEffects & 3%in%selCourseTru){
        
        ccc <- coef(modFin2,which=3,n1=obsPerTra,n2=obsPerTra)$smterms
        
        for(nr in 1:length(levels(dat$g2))){
          
          truth1 <- dat$trueEffHistGame(dat$s,dat$t)*dat$trueEffHistGameFac[nr]
          predEff1 <- ccc[[1]][[nr]]$value
          predEff1[predEff1==0] <- NA
          truth1[lower.tri(truth1)] <- NA
          
          relimseIAGameTru[[nr]] <- sum(c(((predEff1-truth1)^2)),na.rm=T) / sum(c(truth1^2),na.rm = T)
          
        }
        
      }
      
      # Interaction Effect with Random Effect Covariate
      
      
      if(3%in%findEffects && (which(3==findEffects)+1)%in%selCourse){
        
        ccc <- coef(modFin1,which=(which(3==findEffects)+1),n1=obsPerTra,n2=obsPerTra)$smterms
        
        for(nr in 1:length(levels(dat$g3))){
          
          truth1 <- dat$trueEffHistRand[[nr]]
          
          predEff1 <- ccc[[1]][[nr]]$value
          predEff1[predEff1==0] <- NA
          truth1[lower.tri(truth1)] <- NA
          
          relimseIARan[[nr]] <- sum(c(((predEff1-truth1)^2)),na.rm=T)/sum(c(truth1^2),na.rm = T)
          
        }
      }
      
      if(3%in%findEffects && (which(3==findEffects)+1)%in%selCourse){
        
        ccc <- coef(modFin2,which=(which(3==findEffects)+1),n1=obsPerTra,n2=obsPerTra)$smterms
        
        for(nr in 1:length(levels(dat$g3))){
          
          truth1 <- dat$trueEffHistRand[[nr]]
          
          predEff1 <- ccc[[1]][[nr]]$value
          predEff1[predEff1==0] <- NA
          truth1[lower.tri(truth1)] <- NA
          
          relimseIARanTru[[nr]] <- sum(c(((predEff1-truth1)^2)),na.rm=T)/sum(c(truth1^2),na.rm = T)
          
        }
      }
      
      
    }
    
    simDF[[nrSim]] <- (cbind(data.frame(relimseMain=relimseMain, 
                                        relimseIAGame=t(unlist(relimseIAGame)),
                                        relimseIARan=t(unlist(relimseIARan)), 
                                        relimseMainTru=relimseMainTru, 
                                        relimseIAGameTru=t(unlist(relimseIAGameTru)),
                                        relimseIARanTru=t(unlist(relimseIARanTru)), 
                                        mstopIter=mstop(cvr),
                                        nrSim=nrSim),
                             setupDF[i,]))
    
  }
  return(simDF)
  
}, mc.cores = coresSettings)

saveRDS(resSim,file="results/comparison_miss.RDS")
