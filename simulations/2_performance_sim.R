############################################################################
### Simulation code for section 4.2 and 4.6 (step-length comparison)
############################################################################

source("0_libs_funs.R", chdir = T)
if(length(list.files("results")) == 0) dir.create("results")

nrSims = 100
## if you just want to test the code:
if(FALSE) nrSims = 2

### core usage
coresCV = 1
coresSettings = 20


######### settings
addME <- FALSE # TRUE
nuC <- c(0.1#, 1
         )
n <- c(80, 160, 320, 640)
obsPerTra <- c(#20, 
  40#, 60
  )
SNR <- c(1/10, 1, 10)
setup = c("full", "withoutDoubleVar", "histGameIA", "histRandIA")
nrRanEf = 10

######### generate all combinations of different settings
setupDF <- expand.grid(list(setup = setup,
                            nuC = nuC,
                            obsPerTra = obsPerTra,
                            n = n,
                            SNR = SNR,
                            nrRanEf = nrRanEf))

resSim <- vector("list",nrow(setupDF))

# do for all settings
for(i in 1:nrow(setupDF)){
  
  ######### extract settings
  setup = as.character(setupDF$setup[i])
  obsPerTra = setupDF$obsPerTra[i]
  nuC = setupDF$nuC[i]
  n = setupDF$n[i]
  SNR = setupDF$SNR[i]
  nrRanEf = setupDF$nrRanEf[i]
  
  ######### generate data
  dat <- dataGenProc(n = n,
                     obsPerTra = obsPerTra,
                     SNR = SNR,
                     seed = 111,
                     setup = setup,
                     nrOfResp = nrSims,
                     nrRanEf = nrRanEf,
                     nrFacEf = 4
  )
  
  ######### get model specification
  myBlg <-  with(dat,(bols(g2, index=repIDx, df=1) %Xc%
                        brandom(g3, index=repIDx, df=1)))
  
  baseForm <- "Yi ~ 1" 
  terms <- c("bolsc(g2, df=6)",   #1
             "brandomc(g3, df=6)", #2
             "bols(g2, df=2) %Xc% brandom(g3, df=3)",  #3
             "bhistx(X1h, df=15, knots=5, differences=2, standard='length')", #4
             "bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% bolsc(g2, index=repIDx, df=1)",
             "bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% brandomc(g3, index=repIDx, df=1)",
             "bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% myBlg"
  )
  
  ind <- switch(setup,
                histOnly = 4,
                histAndGame = c(4,1),
                histAndRand = c(4,2),
                histGameRand = c(4,1,2,3),
                histGameIA = c(4,5,1),
                histRandIA = c(4,6,2),
                full = c(4,5,6,7,1,2,3),
                withoutDoubleVar = c(4,5,6,1,2,3)
  )
  
  fff <- as.formula(paste(baseForm,paste(terms[ind],collapse = " + "), sep = " + "))
  
  ######### parallelize over repetitions
  simDF <- mclapply(1:nrSims, function(nrSim){
    
    ######### use the ith response matrix
    dat$Yi <- dat$Y[[nrSim]]
    dat$Yvec <- Yvec <- as.vector(dat$Yi)
    
    ######### fit model
    mod2 <- FDboost(fff,
                    timeformula = ~ bbs(t, df=2.5), 
                    data = dat,
                    control = boost_control(mstop = 3500, nu = nuC)
    )
    
    gridEnd = 3500
    gridStart = 1
    
    ######### generate appropriate folds
    ppmat <- createRandomRespFolds(ranVar = dat$g3, sLength = obsPerTra)
    
    ######### validate
    cvr <- if(setup != "histGameIA") cvrisk(mod2, grid = gridStart:gridEnd, 
                  folds = ppmat, mc.cores = coresCV) else
                    cvrisk(mod2, folds = cvLong(id = mod2$id, weights =
                                                  model.weights(mod2), B = 10))
    modFin <- mod2[mstop(cvr)]
    
    findEffects <- which(c(4:7)%in%ind)
    selCourse <- selected(modFin)
    
    relimseMain <- NA
    relimseIAGame <- as.list(rep(NA, length(levels(dat$g2)))) 
    relimseIARan <- as.list(rep(NA, length(levels(dat$g3)))) 
    relimseIAGameRan <- as.list(rep(NA, length(levels(interaction(dat$g2,dat$g3)))))
    
    ######### get relimses
    
    if(1 %in% findEffects & 2 %in% selCourse){
      
      ## Main Effect
      
      trueX1eff <- dat$trueEffHist(dat$s,dat$t)
      trueX1eff[lower.tri(trueX1eff)] <- NA
      predEff <- coef(modFin, which = 2, n1 = obsPerTra, n2 = obsPerTra)$smterms[[1]]$value
      predEff[lower.tri(predEff)] <- NA

      relimseMain <- sum(c(((predEff-trueX1eff)^2)),na.rm = T) / sum(c(trueX1eff^2),na.rm = T)
      
      
    }
    
    if(2 %in% findEffects & 3 %in% selCourse){
      
      ccc <- coef(modFin, which = 3, n1 = obsPerTra, n2 = obsPerTra)$smterms
      
      for(nr in 1:length(levels(dat$g2))){
      
        truth1 <- dat$trueEffHistGame(dat$s,dat$t)*dat$trueEffHistGameFac[nr]
        predEff1 <- ccc[[1]][[nr]]$value
        predEff1[lower.tri(predEff1)] <- NA
        truth1[lower.tri(truth1)] <- NA

        relimseIAGame[[nr]] <- sum(c(((predEff1-truth1)^2)),na.rm = T) / sum(c(truth1^2),na.rm = T)
      
      }
    }
    
    ## Interaction Effect with Random Effect Covariate
    
    
    if(3 %in% findEffects && (which(3 == findEffects) + 1 ) %in% selCourse){
      
      ccc <- coef(modFin, which = (which(3 == findEffects) + 1), n1 = obsPerTra, n2 = obsPerTra)$smterms
      
      for(nr in 1:length(levels(dat$g3))){
        
        truth1 <- dat$trueEffHistRand[[nr]]
        
        predEff1 <- ccc[[1]][[nr]]$value
        predEff1[lower.tri(predEff1)] <- NA
        truth1[lower.tri(truth1)] <- NA
        
        relimseIARan[[nr]] <- sum(c(((predEff1-truth1)^2)),na.rm = T) / sum(c(truth1^2),na.rm = T)
        
      }
    }
    
    ## 3-way Interaction Effect
    
    
    if(4 %in% findEffects & 5 %in% selCourse){
      
      ccc <- coef(modFin, which = 5, n1 = obsPerTra, n2 = obsPerTra)$smterms
      cccSeqFac <- sapply(ccc[[1]][1:(length(ccc[[1]])-2)],
                          function(x)x$add_main)
      cccSeqFac <- gsub("g2=", "", gsub(", g3=", ".", cccSeqFac, fixed=T), fixed=T)
      levIA <- levels(interaction(dat$g2, dat$g3))
      
      for(nr in 1:length(levIA)){
        
        truth1 <- dat$trueDoubleEff[[nr]](dat$s,dat$t)
        
        predEff1 <- ccc[[1]][[which(cccSeqFac == levIA[nr])]]$value
        predEff1[lower.tri(predEff1)] <- NA
        truth1[lower.tri(truth1)] <- NA
        
        relimseIAGameRan[[nr]] <- sum(c(((predEff1-truth1)^2)), na.rm = T) / sum(c(truth1^2), na.rm = T)
        
      }
      
    }
    
    
    return(cbind(data.frame(relimseMain=relimseMain, 
                            relimseIAGame = t(unlist(relimseIAGame)), 
                            relimseIARan = t(unlist(relimseIARan)), 
                            relimseIAGameRan = t(unlist(relimseIAGameRan)),
                            mstopIter = mstop(cvr), nrSim = nrSim), setupDF[i,]))
    
  }, mc.cores = coresSettings)
  
  resSim[[i]] <- simDF
  saveRDS(resSim,file=paste0("results/perf/tempPer",i,".RDS"))
}

######### save results
saveRDS(resSim,file="results/tempPer.RDS")

res <- do.call("rbind",unlist(resSim,recursive=F))

saveRDS(res,file="results/perf_sim.RDS")
