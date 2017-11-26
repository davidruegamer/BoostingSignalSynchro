#########################################################################
#########################################################################
####### Script for reproducing results without 3-way-interaction ########
####### including bootstrap results for aggregated data EEG = Fz ######## 
#######                     and EMG = Fp1-EXG1                   ########
#########################################################################
#########################################################################


# loading libraries
#########################################################################

library(FDboost)
bsplines <- mboost:::bsplines
isMATRIX <- mboost:::isMATRIX

# defining functions
#########################################################################

limitDelta <-  function(s, t) s < t - 3


# loading data
#########################################################################

dAgg <- readRDS("../data/eeg_and_emg_data_aggregated.Rds")


# formatting data
#########################################################################

dAgg$GC_CP <- interaction(dAgg$CP, dAgg$GC)
dAgg$GC_CB <- interaction(dAgg$GC, dAgg$controlBlock)
dAgg$CP_CB <- interaction(dAgg$CP, dAgg$controlBlock)
dAgg$gameCond <- interaction(dAgg$CP, dAgg$GC, dAgg$controlBlock)

dataEx <- list()

dataEx$CP <- matrix(dAgg$CP,ncol=length(unique(dAgg$time)),byrow=T)[,1]
dataEx$GC <- matrix(dAgg$GC,ncol=length(unique(dAgg$time)),byrow=T)[,1]
dataEx$controlBlock <- matrix(dAgg$controlBlock,ncol=length(unique(dAgg$time)),byrow=T)[,1]
dataEx$GC_CP <- matrix(dAgg$GC_CP, ncol=length(unique(dAgg$time)),byrow=T)[,1]
dataEx$GC_CB <- matrix(dAgg$GC_CB, ncol=length(unique(dAgg$time)),byrow=T)[,1]
dataEx$CP_CB <- matrix(dAgg$CP_CB, ncol=length(unique(dAgg$time)),byrow=T)[,1]
dataEx$gameCond <- matrix(dAgg$gameCond, ncol=length(unique(dAgg$time)),byrow=T)[,1]
dataEx$PP <- as.factor(c(matrix(dAgg$vpnr,ncol=length(unique(dAgg$time)),byrow=T)[,1]))
dataEx$PP <- droplevels(dataEx$PP[dataEx$PP!="31"])

nrBoot = 100

subsets <- lapply(1:nrBoot,function(i){ set.seed(i); sample(levels(dataEx$PP), replace=TRUE) })

#########################################################################
# fit all possible models via mclapply
# resps <- c("Fp1_EXG1","EXG2_EXG3","EXG4_EXG5")
# covs <- c("Pz","POz","Fz","FCz")

stillm <- readRDS("results/stillm.RDS")

bootModList <- mclapply(stillm, function(i){
  
  ss <- subsets[[i]]
  
  yVar = "Fp1_EXG1"
  xVar = "Fz"
  
  dataEx$y <- matrix(abs(dAgg[,yVar]),ncol=length(unique(dAgg$time)),byrow=T)
  
  dataEx$X1 <- matrix(dAgg[,xVar],ncol=length(unique(dAgg$time)),byrow=T)
  
  naInd <- is.na(rowMeans(dataEx$X1)) | is.na(rowMeans(dataEx$y))
  indRows <- c(sapply(ss,function(s)which(dataEx$PP==s)))
  
  dataEx$PP <- droplevels(as.factor(dataEx$PP[!naInd][indRows]))
  dataEx$y <- dataEx$y[!naInd,][indRows,]
  dataEx$X1 <- dataEx$X1[!naInd,][indRows,]
  
  dimy <- dim(dataEx$y)[2]
  
  dataEx$CP <- dataEx$CP[!naInd][indRows]
  dataEx$GC <- dataEx$GC[!naInd][indRows]
  dataEx$controlBlock <- dataEx$controlBlock[!naInd][indRows]
  
  dataEx$gameCond <- dataEx$gameCond[!naInd][indRows]
  
  dataEx$GC <- as.factor(dataEx$GC)
  dataEx$CP <- as.factor(dataEx$CP)
  dataEx$controlBlock <- as.factor(dataEx$controlBlock)
  dataEx$GC_CP <- as.factor(dataEx$GC_CP[!naInd][indRows])
  dataEx$GC_CB <- as.factor(dataEx$GC_CB[!naInd][indRows])
  dataEx$CP_CB <- as.factor(dataEx$CP_CB[!naInd][indRows])
  dataEx$gameCond <- as.factor(dataEx$gameCond)
  
  
  
  dataEx$timep <- 1:dim(dataEx$y)[2]
  
  # centering
  dataEx$y <- t(t(dataEx$y) - colMeans(dataEx$y))
  dataEx$X1 <- t(t(dataEx$X1) - colMeans(dataEx$X1))
  
  dataEx$repIDx = rep(1:nrow(dataEx$X1), length(dataEx$timep))
  dataEx$timepVec <- rep(dataEx$timep, each=nrow(dataEx$y))
  dataEx$stimep <- dataEx$timep
  
  condLev <- levels(dataEx$gameCond)
  
  dataExSel <- dataEx
  
  n <- nrow(dataExSel$y)
  nygrid <- ncol(dataExSel$y)
  X1h <- hmatrix(time= rep(dataExSel$timep, each=n), id=rep(1:n, nygrid),
                 x=dataExSel$X1, argvals=dataExSel$stimep,
                 timeLab="timep", idLab="wideIndex", xLab="X1", argvalsLab="stimep")
  
  dataExSel$X1h <- I(X1h)
  

  mod <- FDboost(y ~ 1 + brandomc(PP, df=5) %A% bbs(timep, df=4) +
                   bolsc(CP, df=2.5, intercept=TRUE) %A% bbs(timep, df=8) +
                   bolsc(GC, df=2.5, intercept=TRUE) %A% bbs(timep, df=8) +
                   bolsc(controlBlock, df=2.5, intercept=TRUE)  %A% bbs(timep, df=8)+
                   bolsc(GC_CP, intercept=TRUE, df=2.5) %A% bbs(timep, df=8) +
                   bolsc(GC_CB, intercept=TRUE, df=2.5) %A% bbs(timep, df=8) +
                   bolsc(CP_CB, intercept=TRUE, df=2.5) %A% bbs(timep, df=8) +
                   bolsc(gameCond, intercept=TRUE, df=2.5) %A% bbs(timep, df=8) +
                   bhistx(X1h,
                          limits=limitDelta,
                          df=20, knots=10,
                          differences=2,
                          standard="length"
                   ) +
                   bhistx(X1h,
                          limits=limitDelta,
                          df=5, knots=10,
                          differences=2,
                          standard="length"#,
                          #centerBy=gameCond
                   ) %X%
                   bolsc(gameCond, df=4, intercept=TRUE, index=repIDx) +
                   bhistx(X1h,
                          limits=limitDelta,
                          df=5, knots=10,
                          differences=2,
                          standard="length"#,
                          #centerBy=PP
                   ) %X%
                   brandomc(PP, df=4, index=repIDx),
                 control=boost_control(mstop=5000, trace=TRUE),
                 timeformula = ~ bbs(timep),
                 data = dataExSel
  )

  saveRDS(mod,file=paste0("v17_apr16_agg_mod_allpp_mstop2500_knots10_wo31_Fz_Fp1EXG1_nr",i,".Rds"))

  # get the optimal stopping iteration per model
  
  lld <- length(levels(dataEx$PP))
  fff = 6
  
  ppVec <- dataEx$PP
  ppVec2 <- rep(dataEx$PP, each=384)#ncol(dataExSel$y))
  nrFolds <- fff
  nrPPperFold <- floor(lld/nrFolds)
  
  ppMat2 <- matrix(rep(ppVec2,nrFolds),ncol=nrFolds)
  ppIn <- split(levels(dataEx$PP),1:fff)
  ppMat2 <- sapply(0:(nrFolds-1),function(i)
    as.numeric(sapply(ppMat2[,i+1],function(x)!x%in%ppIn[[i+1]])))
  
  try(cvr <- cvrisk(mod,grid=1:5000,
                    folds=ppMat2,
                    mc.cores=1))
  try(print(paste(xVar,yVar,fff,mstop(cvr),sep="->")))
  saveRDS(cvr,
          paste0("v17_apr16_cv_all_10_wo31_Fz_Fp1EXG1_6_fold_nr",i,".Rds"))
  
},mc.cores=5)

