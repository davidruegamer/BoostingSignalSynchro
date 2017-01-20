############################################################################
### Function to simulate data for a main historical effect, random or
### factor-specific historical effects as well as double-varying hist. eff.
### Code partly adopted from simulation by Sarah Brockhaus
############################################################################

dataGenProc <- function(n = 320, # number of trajectories 
                        obsPerTra = 40, # number of observations
                        SNR = 10, # signal-to-noise ratio
                        seed = NULL, # seed
                        seedRanEf = 1337, # seed for random effect
                        seedRanEfHist = 1337, # seed for random historical effect generation
                        seedFacEf = 1337, # seed for factor effects
                        setup = c("histOnly",    # how data is generated
                                  "histAndGame",
                                  "histAndRand",
                                  "histGameRand",
                                  "histGameIA",
                                  "histRandIA",
                                  "full",
                                  "withoutDoubleVar"),
                        k = 15, # df for spline on which covariates are generated
                        meFun = function(t)(sin(1500*t)/20), # measurement error function
                        addME = FALSE, # whether to add additional oscillating measurement error
                        nrOfResp = 1, # how many data sets should be generated (same covariates, different response)
                        nrRanEf = 8, # number of random effect groups > 1 (no impact for settings w/o random effect)
                        # partRanEf = 16, # defines the number of observations per random effect
                        nrFacEf = 2, # number of factor groups
                        redFacForRanEf  = 10, # random coefficient surfaces are divided by this number
                        intFun = function(t) cos(t),
                        trueEffHist = function(s,t){ # true historical effect
                          
                          ret <- outer(s, t, function(s, t) 2*sin(sqrt(s*t+0.5)*pi))
                          
                          for(i in 1:length(s)){
                            for(j in 1:length(t)){
                              if(s[i] > t[j]){
                                ret[i, j] <- 0
                              } 
                            }
                          }
                          return(ret)
                        },
                        trueEffHistGame = function(s,t){ # true game specific historical effect
                          
                          ret <- outer(s, t, function(s, t) s*cos(pi*sqrt(t)) )
                          
                          for(i in 1:length(s)){
                            for(j in 1:length(t)){
                              if(s[i] > t[j]){
                                ret[i, j] <- 0
                              } 
                            }
                          }
                          return(ret)
                        }
)
{
  
  if(!is.null(seed)) set.seed(seed)
  if(n%%nrFacEf!=0 | n%%nrRanEf!=0) # implementation for an unbalanced design not possible in a generic manner
    stop("Please specify n as a multiple of nrRanEf and nrFacEf. Otherwise, the setting is unbalanced (not supported).")
  
  ################## create the basic objects #####################
  
  setup <- match.arg(setup)
  x1 <- runif(n)
  x1orig <- x1
  x1 <- x1 - mean(x1)
  x2 <- runif(n)
  x2 <- x2 - mean(x2)
  
  
  
  g2 <- gl(n=nrFacEf, k=n/nrFacEf, length=n)
  g3 <- factor(rep(1:nrRanEf,n/nrRanEf))
  
  rf <- function(s=seq(0,1,length=100), k=k, center=FALSE) {
    drop(ns(s, k, int=TRUE) %*% runif(k, -3, 3)+15)
  }
  s <- seq(0, 1, length = obsPerTra)
  t <- seq(0, 1, length = obsPerTra)
  X1 <- t(replicate(n, rf(s, k=k)))
  X1cen <- t(t(X1) - colMeans(X1))

  ################## evaluate the given coefficient functions #####################
  
  L <- integrationWeightsLeft(X1=X1cen, xind=s)  # Riemann integration weights      
  
  mainBeta <- function(X1, s, t){
    (L*X1) %*% ( trueEffHist(s,t))
  }
  g2Beta <- function(X1, s, t){
    (L*X1) %*% ( trueEffHistGame(s,t))
  }
  g2g3BetaVar <- function(X1, s, t, betaFac){ # for random factor-specific hist effects
    (L*X1) %*% ( betaFac*trueEffHistGame(s,t))
  }
  
  ################## evaluate all coefficient functions #####################
  
  funxbetaMain <-  mainBeta(X1cen, s, t)
  funxbetaG2 <-  g2Beta(X1cen, s, t)

  # random part
  
  set.seed(seedFacEf)
  
  facEff <- scale(runif(nrFacEf,-5,5))
  
  set.seed(seedRanEf)
  
  ranSim <- c(scale(rnorm(nrRanEf)))
  ranSim2 <- t(facEff%*%ranSim)
  
  raneffs <- lapply(genRanHist(s, t, nrLevels = nrRanEf, seed = seedRanEfHist),function(r)r/redFacForRanEf)
  funxRanHist <- lapply(raneffs,function(beta)(L*X1cen)%*%beta)

  funxRanGameHist <- lapply(1:(nrRanEf*nrFacEf), function(i) g2g3BetaVar(X1=X1cen, s=s, t=t, betaFac=c(t(ranSim2[i]))))
  
  ################## create the true underlying response matrix #####################
  
  Ytrue <-  
    matrix(rep(intFun(t),nrow(X1cen)),ncol=length(t),byrow=T) + # intercept
    
    as.numeric(!setup%in%c("histOnly","histAndRand","histRandIA"))* # game condition effect
    (
      Reduce("+",lapply(1:nrFacEf,function(i)
        matrix(rep(facEff[i]*t,nrow(X1cen)),ncol=length(t),byrow=T)*as.numeric(g2==i)
      ))
    ) + 
    
    as.numeric(!setup%in%c("histOnly","histAndGame","histGameIA"))*( # random effect
      Reduce("+",lapply(1:nrRanEf,function(i)
        matrix(rep(ranSim[i]*t,nrow(X1cen)),ncol=length(t),byrow=T)*as.numeric(g3==i)
      ))
    ) + 
    
    as.numeric(setup%in%c("histGameRand","full","withoutDoubleVar"))*( # random interaction effect
      Reduce("+",lapply(1:(nrRanEf*nrFacEf),
                        function(i)matrix(rep(c(t(ranSim2))[i]*t,nrow(X1cen)),ncol=length(t),byrow=T)*
                          as.numeric(interaction(g2,g3)==levels(interaction(g2,g3))[i])))
    ) +
    
    funxbetaMain + # main  historical effect
    
    as.numeric(!setup%in%c("histOnly","histAndRand","histRandIA"))* # game condition specific hist
    (
      Reduce("+",lapply(1:nrFacEf,function(i) (funxbetaG2 * facEff[i]) * as.numeric(g2==i)))
    ) +
    
    as.numeric(!setup%in%c("histOnly","histAndGame","histGameIA","histAndRand"))*  # random historical effect
    (Reduce("+",lapply(1:nrRanEf,function(i)funxRanHist[[i]] * as.numeric(g3==i)))) + 
    
    as.numeric(setup=="full")*( # double varying historical effect
      Reduce("+", lapply(1:(nrRanEf*nrFacEf), function(i) 
        funxRanGameHist[[i]] * 
          as.numeric(interaction(g2,g3)==levels(interaction(g2,g3))[i])))
    )
  
  ################## create list of all components #####################
  
  dat <- list(x1 = x1, x1orig = x1orig, x2 = x2, g2 = g2, 
              g3 = g3, X1 = I(X1), X1cen = I(X1cen), s = s, t = t)
  dat$trueEffHist <- trueEffHist
  dat$trueEffHistGame <- trueEffHistGame
  dat$trueEffHistGameFac <- facEff
  dat$trueEffHistRand <- raneffs
  dat$trueDoubleEff <- lapply(1:(nrRanEf*nrFacEf),function(i)function(s,t)c(t(ranSim2[i]))*trueEffHistGame(s,t))
  
  # objects needed for vector representations in FDboost
  dat$repIDx <- repIDx <- rep(1:nrow(dat$X1), length(dat$t))
  dat$tvec <- tvec <- rep(t, each=nrow(dat$X1)) 
  
  ################## hmatrix structure for historical effect #####################
  
  X1h <- hmatrix(time = tvec, id = repIDx, 
                 x = X1cen, argvals = s, 
                 timeLab="t", idLab="wideIndex", xLab="X1h", argvalsLab="s")
  
  dat$X1h <- I(X1h)   
  
  ################## create response variables #####################
  
  if(addME) Ytrue <- Ytrue + meFun(t) # measurement error
  dat$Ytrue <- Ytrue
  
  if(SNR==Inf) dat$Y <- dat$Yvec <- Ytrue
  
  if(nrOfResp>1){ # generate multiple response matrices
    
    Y <- vector("list", nrOfResp)
    Y <- lapply(1:nrOfResp,function(i)
      Ytrue + sd(as.vector(Ytrue))/sqrt(SNR) * matrix(scale(rnorm(n * obsPerTra)), nrow = n)
    )
    dat$Y <- Y
    
  }else{ # generate only one response matrix
    
    Y <- Ytrue + sd(as.vector(Ytrue))/sqrt(SNR) * matrix(scale(rnorm(n * obsPerTra)), nrow = n)
    
    dat$Y <- Y
    dat$Yvec <- Yvec <- as.vector(dat$Y)
    
  }
  
  return(dat)
  
}