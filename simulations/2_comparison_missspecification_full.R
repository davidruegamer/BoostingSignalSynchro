############################################################################
### Simulation code for section 4.2 (missspecification)
############################################################################

library("mboost")
library("FDboost")
library("splines")        
source("../functions/ranHist.R")
source("../functions/createRandomRespFolds.R")
source("../functions/dataGenProc.R")
source("../functions/matggplot.R")

######### settings
addME <- TRUE
nuC <- c(0.1)
n <- c(80, 320)
obsPerTra <- c(40) #, 60)
SNR <- c(1/20, 5, 20)
setup = c("full", "withoutDoubleVar")
nrSims = 100

######### generate all combinations of different settings
setupDF <- expand.grid(list(setup=setup,
                            n=n,
                            SNR=SNR))
setupDF$setup <- as.character(setupDF$setup)

######### parallelize over different settings
resSim <- mclapply(1:nrow(setupDF),function(i){
  
  ######### extract settings
  setup = setupDF$setup[i]
  n = setupDF$n[i]
  SNR = setupDF$SNR[i]
  
  ######### generate data
  dat <- dataGenProc(n = n,
                     obsPerTra = obsPerTra,
                     SNR = SNR,
                     seed = 12,
                     setup = "full",
                     nrOfResp = nrSims,
                     nrRanEf = 8
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
  
  # fff <- as.formula(paste(baseForm,paste(terms[ind],collapse=" + "),sep=" + "))
  
  simDF <- vector("list",nrSims)
  
  ######### do for nrSims repetitions
  for(nrSim in 1:nrSims){
    
    obsPerTra <- c(40)
    
    ######### use the ith response matrix
    dat$Yi <- if(nrSims>1) dat$Y[[nrSim]] else dat$Y
    dat$Yvec <- Yvec <- as.vector(dat$Yi)
    
    ######### fit model
    myBlg <<-  with(dat,(bols(g2, index=repIDx, df=1) %Xc%
                           brandom(g3, index=repIDx, df=1)))
    
    
    if(setup=="full"){
      
      mod2 <- FDboost(Yi ~ 1 +  
                        bhistx(X1h, df=15, knots=5, differences=2, standard='length') +
                        bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% bolsc(g2, index=repIDx, df=1) +
                        bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% brandomc(g3, index=repIDx, df=1) +
                        bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% myBlg +
                        bolsc(g2, df=6) +
                        brandomc(g3, df=6) +
                        bols(g2, df=2) %Xc% brandom(g3, df=3),
                      timeformula = ~ bbs(t, df=2.5), data=dat,
                      control=boost_control(mstop=2500,nu=0.1)
      )
      
    }else{
      
      mod2 <- FDboost(Yi ~ 1 +
                        bhistx(X1h, df=15, knots=5, differences=2, standard='length') +
                        bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% bolsc(g2, index=repIDx, df=1) +
                        bhistx(X1h, df=15, knots=5, differences=2, standard='length') %X% brandomc(g3, index=repIDx, df=1) +
                        bolsc(g2, df=6) +
                        brandomc(g3, df=6) +
                        bols(g2, df=2) %Xc% brandom(g3, df=3) ,
                      timeformula = ~ bbs(t, df=2.5), data=dat,
                      control=boost_control(mstop=2500,nu=0.1)
      )
      
    }
    
    gridEnd = 2500
    gridStart = 1
    
    ######### generate appropriate folds
    ppmat <- createRandomRespFolds(ranVar = dat$g3, sLength = obsPerTra)
    
    ######### validate
    cvr <- cvrisk(mod2, grid=gridStart:gridEnd, 
                  folds=ppmat, mc.cores=1)#4)
    modFin <- mod2[mstop(cvr)]
    
    findEffects <- which(c(4:7)%in%ind)
    selCourse <- selected(modFin)
    
    relimseMain <- relimseIAGame <- NA
    relimseIARan <- as.list(rep(NA,length(levels(dat$g3)))) 
    # relimseIAGameRan <- as.list(rep(NA,length(levels(interaction(dat$g2,dat$g3)))))
    
    ######### get relimses
    
    
    if(1%in%findEffects & 2%in%selCourse){
      
      ## Main Effect
      
      trueX1eff <- dat$trueEffHist(dat$s,dat$t)
      trueX1eff[trueX1eff==0] <- NA
      predEff <- coef(modFin,which=2,n1=obsPerTra,n2=obsPerTra)$smterms[[1]]$value
      predEff[predEff==0] <- NA
      
      # print(matggplot(listOfMatrices = list(trueX1eff,predEff), 
      #                 namesForGrids = c("truth","FDboost")))
      # 
      relimseMain <- sum(c(((predEff-trueX1eff)^2)),na.rm=T)/sum(c(trueX1eff^2),na.rm = T)
      
      
    }
    
    if(2%in%findEffects & 3%in%selCourse){
      
      ccc <- coef(modFin,which=3,n1=obsPerTra,n2=obsPerTra)$smterms
      
      truth1 <- dat$trueEffHistGame(dat$s,dat$t)
      predEff1 <- ccc[[1]][[1]]$value
      predEff1[predEff1==0] <- NA
      truth1[lower.tri(truth1)] <- NA
      #truth1[sapply(c(truth1),function(s)all.equal(s,0)=="TRUE")] <- NA
      
      print(matggplot(listOfMatrices = list(truth1,predEff1),
                      namesForGrids = c("truth","FDboost")))
      
      relimseIAGame <- sum(c(((predEff1-truth1)^2)),na.rm=T)/sum(c(truth1^2),na.rm = T)
      
    }
    
    ## Interaction Effect with Random Effect Covariate
    
    
    if(3%in%findEffects && (which(3==findEffects)+1)%in%selCourse){
      
      ccc <- coef(modFin,which=(which(3==findEffects)+1),n1=obsPerTra,n2=obsPerTra)$smterms
      
      for(nr in 1:length(levels(dat$g3))){
        
        truth1 <- dat$trueEffHistRand[[nr]]
        
        predEff1 <- ccc[[1]][[nr]]$value
        predEff1[predEff1==0] <- NA
        truth1[lower.tri(truth1)] <- NA
        #truth1[sapply(c(truth1),function(s)all.equal(s,0)=="TRUE")] <- NA
        
        # print(matggplot(listOfMatrices = list(truth1,predEff1), 
        #                 namesForGrids = c("truth","FDboost")))
        
        relimseIARan[[nr]] <- sum(c(((predEff1-truth1)^2)),na.rm=T)/sum(c(truth1^2),na.rm = T)
        
      }
    }
    
    # ## 3-way Interaction Effect
    # 
    # 
    # if(4%in%findEffects & 5%in%selCourse){
    #   
    #   ccc <- coef(modFin,which=5,n1=obsPerTra,n2=obsPerTra)$smterms
    #   cccSeqFac <- sapply(ccc[[1]][1:(length(ccc[[1]])-2)],
    #                       function(x)x$add_main)
    #   cccSeqFac <- gsub("g2=","",gsub(", g3=",".",cccSeqFac,fixed=T),fixed=T)
    #   levIA <- levels(interaction(g2,g3))
    #   
    #   for(nr in 1:length(levIA)){
    #     
    #     truth1 <- dat$trueDoubleEff[[nr]](dat$s,dat$t)
    #     
    #     predEff1 <- ccc[[1]][[which(cccSeqFac==levIA[nr])]]$value
    #     predEff1[predEff1==0] <- NA
    #     truth1[lower.tri(truth1)] <- NA
    #     
    #     print(matggplot(listOfMatrices = list(truth1,predEff1), 
    #                     namesForGrids = c("truth","FDboost")))
    #     
    #     relimseIAGameRan[[nr]] <- sum(c(((predEff1-truth1)^2)),na.rm=T)/sum(c(truth1^2),na.rm = T)
    #     
    #   }
    #   
    # }
    
    
    simDF[[nrSim]] <- (cbind(data.frame(relimseMain=relimseMain, relimseIAGame=relimseIAGame, 
                                        relimseIARan=t(unlist(relimseIARan)), 
                                        #relimseIAGameRan=t(unlist(relimseIAGameRan)),
                                        mstopIter=mstop(cvr),nrSim=nrSim),setupDF[i,]))
    
  }
  return(simDF)
  
},mc.cores=12)

saveRDS(resSim,file="tempMis.RDS")

res <- do.call("rbind",unlist(resSim,recursive=F))

saveRDS(res,file="results/misspecification_sim.RDS")
