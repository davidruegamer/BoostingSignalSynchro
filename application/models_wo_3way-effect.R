#########################################################################
#########################################################################
####### Script for reproducing results without 3-way-interaction ########
#######                   and aggregated data                    ########
#########################################################################
#########################################################################


# loading libraries
#########################################################################

library(FDboost)
source("../functions/createRandomRespFolds.R")
library(ggplot2)
library(RColorBrewer)
palSpec <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
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



# fit all possible models via mclapply
#########################################################################

resps <- c("Fp1_EXG1","EXG2_EXG3","EXG4_EXG5")
covs <- c("Pz","POz","Fz","FCz")

predictionList <- mclapply(as.list(resps),function(yVar){

  predListXvar <- mclapply(as.list(covs), function(xVar){

# yVar = "Fp1_EXG1"
# xVar = "Fz"

    dataEx$y <- matrix(abs(dAgg[,yVar]),ncol=length(unique(dAgg$time)),byrow=T)

    dataEx$X1 <- matrix(dAgg[,xVar],ncol=length(unique(dAgg$time)),byrow=T)

    naInd <- is.na(rowMeans(dataEx$X1)) | is.na(rowMeans(dataEx$y)) | dataEx$PP%in%c("31")#,"27","28","33") 

    if(sum(naInd)!=0){

      dataEx$PP <- droplevels(dataEx$PP[!naInd])
      dataEx$y <- dataEx$y[!naInd,]
      dataEx$X1 <- dataEx$X1[!naInd,]

      dimy <- dim(dataEx$y)[2]

      dataEx$CP <- dataEx$CP[!naInd]
      dataEx$GC <- dataEx$GC[!naInd]
      dataEx$controlBlock <- dataEx$controlBlock[!naInd]

      dataEx$gameCond <- dataEx$gameCond[!naInd]

      dataEx$GC <- as.factor(dataEx$GC)
      dataEx$CP <- as.factor(dataEx$CP)
      dataEx$controlBlock <- as.factor(dataEx$controlBlock)
      dataEx$GC_CP <- as.factor(dataEx$GC_CP[!naInd])
      dataEx$GC_CB <- as.factor(dataEx$GC_CB[!naInd])
      dataEx$CP_CB <- as.factor(dataEx$CP_CB[!naInd])

    }

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
                     #                  #                           bols(gameCond, df=4,
                     #                  #                                # contrasts.arg = "contr.dummy",
                     #                  #                                index=repIDx
                     #                  #                           ) + # %X%
                     brandomc(PP, df=4, index=repIDx),
                   control=boost_control(mstop=5000, trace=TRUE),
                   timeformula = ~ bbs(timep),
                   data = dataExSel
    )
    
    saveRDS(mod,file=paste0("v17_feb16_agg_mod_allpp_mstop5000_knots10_wo31_",xVar,"_",yVar,".Rds"))
    
  },mc.cores=4)
},mc.cores=3)


# get the optimal stopping iteration per model
#########################################################################
# 

dataEx$PP <- as.factor(c(matrix(dAgg$vpnr,ncol=length(unique(dAgg$time)),byrow=T)[,1]))
dataEx$PP <- droplevels(dataEx$PP[dataEx$PP!="31"])

for(fff in c(6,12)){

  ppVec <- dataEx$PP
  ppVec2 <- rep(dataEx$PP, each=384)#ncol(dataExSel$y))
  nrFolds <- fff
  nrPPperFold <- round(length(levels(dataEx$PP))/nrFolds)

ppMat2 <- matrix(rep(ppVec2,nrFolds),ncol=nrFolds)
ppMat2 <- createRandomRespFolds(ranVar = ppVec, folds = fff, sLength = 384)

mclapply(as.list(resps),function(yVar){
  mclapply(as.list(covs), function(xVar){


    mod <- readRDS(paste0("v17_feb16_agg_mod_allpp_mstop5000_knots10_wo31_",xVar,"_",yVar,".Rds"))

    try(cvr <- cvrisk(mod,grid=1:5000,
                           folds=ppMat2,
                           mc.cores=12))
    try(print(paste(xVar,yVar,fff,mstop(cvr),sep="->")))
    saveRDS(cvr,
            paste0("v17_feb16_cv_all_10",xVar,"_",yVar,"_",fff,
                   "folds_corrected_grid5000_wo31.Rds"))

  },mc.cores=1)
},mc.cores=3)

}


for(yVar in resps){

  for(xVar in covs){

    mod <- readRDS(paste0("v17_feb16_agg_mod_allpp_mstop5000_knots10_wo31_",xVar,"_",yVar,".Rds"))

    for(fff in c(6,12)){

      mss <- mstop(readRDS(paste0("v17_feb16_cv_all_10",xVar,"_",yVar,"_",fff,
                                  "folds_corrected_grid5000_wo31.Rds")))
      modMs <- mod[mss]
      pred <- coef(mod, which=11)
      predMain <- coef(mod, which=10)$smterms[[1]]$value
      preds <- pred$smterms[[1]][1:(pred$smterms[[1]]$numberLevels)]
      predX <- lapply(preds,"[[","x")
      predY <- lapply(preds,"[[","y")
      predVals <- lapply(lapply(preds,"[[","value"),function(x)x+predMain)
      namesPreds <- levels(sapply(preds,"[[","z"))

      predVals <- lapply(predVals,function(mat){mat[mat==0]<-NA;return(mat)})

      res <- data.frame(z=c(sapply(predVals,function(z)c((z)))),
                        x=c(sapply(predX,function(x)rep(x,each=length(x)))),
                        y=c(sapply(predY,function(y)rep(y,times=length(y)))),
                        cond=rep(namesPreds,each=length(predX[[1]])^2))

      g <- ggplot(res, aes(x=x,y=y,z=z,fill=z)) + geom_tile() + theme_bw() + facet_grid(cond~.)

      g <- g + scale_fill_gradientn(colours = palSpec(250)) + stat_contour(binwidth = 2) +
        ggtitle(paste0(yVar," explained by ",xVar," with ",fff," folds"))

      print(g)

      saveRDS(g, file=paste0("addedEff_v17_feb16_cv_all_10",
                             xVar,"_",yVar,"_",fff,"folds_corrected_grid5000_wo31.Rds"))

    }
  }
}
