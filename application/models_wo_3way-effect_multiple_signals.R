#########################################################################
#########################################################################
#######    Script for fitting a model on all EEG-signals with    ########
#######                     aggregated data                      ########
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

dAgg1 <- readRDS(file="../data/eeg_and_emg_data_all_cov_aggregated_part1.Rds")
dAgg2 <- readRDS(file="../data/eeg_and_emg_data_all_cov_aggregated_part2.Rds")
dAgg <- cbind(dAgg1, dAgg2[,c(-1*c(1:8))])

# formatting data
#########################################################################

dAgg$GC_CP <- interaction(dAgg$CP, dAgg$GC)
dAgg$GC_CB <- interaction(dAgg$GC, dAgg$controlBlock)
dAgg$CP_CB <- interaction(dAgg$CP, dAgg$controlBlock)
dAgg$gameCond <- interaction(dAgg$CP, dAgg$GC, dAgg$controlBlock)

dAgg <- dAgg[order(dAgg$vpnr,dAgg$CP,dAgg$GC,dAgg$controlBlock,dAgg$time),]

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

yVar = "Fp1_EXG1"

dataEx$y <- matrix(abs(dAgg[,yVar]),ncol=length(unique(dAgg$time)),byrow=T)

covarMats <- lapply(colnames(dAgg)[-1*c(1:8,73:76)], function(xVar)
  matrix(dAgg[,xVar],ncol=length(unique(dAgg$time)),byrow=T))

names(covarMats) <- paste0("X_", colnames(dAgg)[-1*c(1:8,73:76)])

dataEx <- append(dataEx, covarMats)

naInd <- dataEx$PP%in%c("31")

if(sum(naInd)!=0){
  
  dataEx$PP <- droplevels(dataEx$PP[!naInd])
  dataEx$y <- dataEx$y[!naInd,]
  for(nn in names(covarMats)){
    
    dataEx[[nn]] <- dataEx[[nn]][!naInd,]
    
  }
  
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
for(nn in names(covarMats)){
  
  dataEx[[nn]] <- t(t(dataEx[[nn]]) - colMeans(dataEx[[nn]]))
  
}

dataEx$repIDx = rep(1:nrow(dataEx$y), length(dataEx$timep))
dataEx$timepVec <- rep(dataEx$timep, each=nrow(dataEx$y))
dataEx$stimep <- dataEx$timep

condLev <- levels(dataEx$gameCond)

dataExSel <- dataEx

n <- nrow(dataExSel$y)
nygrid <- ncol(dataExSel$y)

for(nn in names(covarMats)){
  
  X1h <- hmatrix(time= rep(dataExSel$timep, each=n), id=rep(1:n, nygrid),
                 x=dataExSel[[nn]], argvals=dataExSel$stimep,
                 timeLab="timep", idLab="wideIndex", xLab=paste0("X_",nn), argvalsLab="stimep")
  
  dataExSel[[paste0("Xh_",nn)]] <- I(X1h)

}

mod <- FDboost(
  as.formula(paste0(
  "y ~ 1 + brandomc(PP, df=5) %A% bbs(timep, df=4) +
                 bolsc(CP, df=2.5, intercept=TRUE) %A% bbs(timep, df=8) +
                 bolsc(GC, df=2.5, intercept=TRUE) %A% bbs(timep, df=8) +
                 bolsc(controlBlock, df=2.5, intercept=TRUE)  %A% bbs(timep, df=8)+
                 bolsc(GC_CP, intercept=TRUE, df=2.5) %A% bbs(timep, df=8) +
                 bolsc(GC_CB, intercept=TRUE, df=2.5) %A% bbs(timep, df=8) +
                 bolsc(CP_CB, intercept=TRUE, df=2.5) %A% bbs(timep, df=8) +
                 bolsc(gameCond, intercept=TRUE, df=2.5) %A% bbs(timep, df=8) +",
  paste(paste0("bhistx(Xh_", names(covarMats),
                        ", limits=limitDelta,
                        df=20, knots=10,
                        differences=2,
                        standard='length'
                 ) +
                 bhistx(Xh_", names(covarMats),
                        ", limits=limitDelta,
                        df=5, knots=10,
                        differences=2,
                        standard='length'
                 ) %X%
                 bolsc(gameCond, df=4, intercept=TRUE, index=repIDx) +
                 bhistx(Xh_",names(covarMats),
                        ", limits=limitDelta,
                        df=5, knots=10,
                        differences=2,
                        standard='length'
                 ) %X%
                 brandomc(PP, df=4, index=repIDx)"), collapse = " + "))),
               control=boost_control(mstop=2000, trace=TRUE),
               timeformula = ~ bbs(timep),
               data = dataExSel
)

saveRDS(mod,file=paste0("allCovarMod.Rds"))
mod <- readRDS("~/Downloads/BoostingSignalSynchro-master/application/allCovarMod.Rds")
saveRDS(table(selected(mod)), file="selFreqAllCovarMod.RDS")

resSqMat <- matrix(mod$resid()^2, nrow=184)
yMat <- matrix(mod$response, nrow=184)

Rsq <- pmax(0, 1 - (colSums(resSqMat) / colSums((yMat - colMeans(yMat))^2)))
saveRDS(Rsq, file="allCovarRsq.RDS")
