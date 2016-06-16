############################################################################
### Simulation code for section 4.2 and 4.4 (step-length comparison)
############################################################################

source("0_libs_funs.R", chdir = T)
if(length(list.files("results")) == 0) dir.create("results")

nrSims = 100
## if you just want to test the code:
if(FALSE) nrSims = 2

### core usage
coresCV = 5
coresSettings = 5


######### settings
addME <- FALSE # TRUE
nuC <- c(0.1)
n <- c(160) #, 200)
obsPerTra <- c(60, 180, 380)
SNR <- c(1/100, 1/10, 1)
setup = c("histGameIA")

######### generate all combinations of different settings
setupDF <- expand.grid(list(#setup=setup,
                            obsPerTra=obsPerTra,
                            #n=n,
                            SNR=SNR))

resSim <- vector("list",nrow(setupDF))

# do for all settings
for(i in 1:nrow(setupDF)){
  
  ######### extract settings
  #setup = as.character(setupDF$setup[i])
  obsPerTra = setupDF$obsPerTra[i]
  #nuC = setupDF$nuC[i]
  #n = setupDF$n[i]
  SNR = setupDF$SNR[i]
  
  ######### generate data
  dat <- dataGenProc(n = n,
                     obsPerTra = obsPerTra,
                     SNR = SNR,
                     seed = 111,
                     setup = setup,
                     nrOfResp = nrSims,
                     nrFacEf = 8
  )

 ######### parallelize over repetitions
  simDF <- mclapply(1:nrSims,function(nrSim){
    
    ######### use the ith response matrix
    dat$Yi <- dat$Y[[nrSim]]
    dat$Yvec <- Yvec <- as.vector(dat$Yi)
    
    ######### fit model
    mod2 <- FDboost(Yi ~ 1 + bhistx(X1h, df=15, knots=5, differences=2, standard='length') + 
                      bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% 
                      bolsc(g2, index=repIDx, df=1) + 
                      bolsc(g2, df=6),
                    timeformula = ~ bbs(t, df=2.5), data=dat,
                    control=boost_control(mstop=1500,nu=nuC)
    )
    
    gridEnd = 1500
    gridStart = 1
    
    ######### validate
    cvr <- cvrisk(mod2, grid=gridStart:gridEnd, mc.cores = coresCV, 
                  folds = cvLong(id = mod2$id, weights = model.weights(mod2), B = 10))
    modFin <- mod2[mstop(cvr)]
    
    findEffects <- which(c(4:7)%in%c(4,5,1))
    selCourse <- selected(modFin)
    
    relimseMain <- relimseIAGame <- NA

    ######### get relimses
    
    
    if(1%in%findEffects & 2%in%selCourse){
      
      ## Main Effect
      
      trueX1eff <- dat$trueEffHist(seq(0,1,l=40),seq(0,1,l=40))
      trueX1eff[trueX1eff==0] <- NA
      predEff <- coef(modFin,which=2,n1=40,n2=40)$smterms[[1]]$value
      predEff[predEff==0] <- NA
      
      relimseMain <- sum(c(((predEff-trueX1eff)^2)),na.rm = T)/sum(c(trueX1eff^2),na.rm = T)
      
      
    }
    
    if(2%in%findEffects & 3%in%selCourse){
      
      ccc <- coef(modFin,which=3,n1=40,n2=40)$smterms
     
      relimseIAGame <- median(sapply(1:8,function(j){ 
        
        truth1 <- dat$trueEffHistGame(seq(0,1,l=40),seq(0,1,l=40))*dat$trueEffHistGameFac[j]
        predEff1 <- ccc[[1]][[j]]$value
        predEff1[predEff1==0] <- NA
        truth1[lower.tri(truth1)] <- NA
        
        sum(c(((predEff1-truth1)^2)),na.rm = T)/sum(c(truth1^2),na.rm = T)
        
      }))
      
    }
    
    return(cbind(data.frame(relimseMain=relimseMain, relimseIAGame=relimseIAGame, 
                                        mstopIter=mstop(cvr),nrSim=nrSim),setupDF[i,]))
    
  }, mc.cores = coresSettings)
  
  resSim[[i]] <- simDF
  
}

######### save results
saveRDS(resSim,file="results/tempPer_obsPerTra.RDS")

res <- do.call("rbind",unlist(resSim,recursive=F))

saveRDS(res,file="results/perf_sim_obsPerTra.RDS")
