############################################################################
### Simulation code for section 4.1 
############################################################################

source("0_libs_funs.R", chdir = T)
if(length(list.files("results")) == 0) dir.create("results")

nrSims = 100
## if you just want to test the code:
if(FALSE) nrSims = 2

### core usage
coresCV = 5
coresSettings = 9

# true underlying function for simulation of multi-modal coefficient surfaces
hillyFun <- function(s,t) sin(10 * abs(s-t)) * cos(10 * t)
# limit function
lD <- function(s, t) s < t
# fixed settings
obsPerTra = 40


nn <- c(80, 160, 320, 640)
ss <- c(1/10, 1, 10)
bd <- c(#5, 7, 9, 
  11) 
k <- c(#5, 
  15)
setupDF <- expand.grid(list(nn = nn, ss = ss, bd = bd, k = k))

set.seed(42)

# parallelize simulations over different settings

res <- mclapply(1:nrow(setupDF), function(set){
  
  # extract settings
  nn = setupDF[set,]$nn
  ss = setupDF[set,]$ss
  bd = setupDF[set,]$bd
  k = setupDF[set,]$k
  
  resDF <- vector("list", nrSims)
  
  # for nrSims repetitions do the following
  for(nrSim in 1:nrSims){
    
    # 84, 21, 52, 53
    
    obsPerTra = 40
    
    # create data using a variation of the pffrSim function
    data1 <- pffrSimVar(scenario=c("int", "ff"), 
                        n = nn, 
                        nxgrid = obsPerTra, 
                        nygrid = obsPerTra, 
                        limits = lD, 
                        SNR = ss,
                        simFun = hillyFun, 
                        bd = bd)
    
    # center data
    data1$X1 <- t(t(data1$X1) - colMeans(data1$X1))
    
    X1 <- matrix(data1$X1, ncol = nn)
    Y <- matrix(data1$Y, ncol = nn)
    
    test1 <- hillyFun
    t <- attr(data1, "yindex")
    s <- attr(data1, "xindex")
    
    # calcuate true coefficient function
    trueBeta <- outer(s, t, test1) * outer(s, t, lD)
    trueBeta[trueBeta == 0] <- NA
    
    # fit pffr function
    m1 <- pffr(Y ~ 1 + ff(X1, xind=s, 
                          splinepars = list(bs = "ps", m = list(2,1), k = k), 
                          yind = t,
                          limits = lD), 
               data = data1)
    p <- plot(m1, pers = FALSE, select = 2, 
              scheme = 2, rug = FALSE, n2 = obsPerTra)

    data1 <- as.list(data1)
    data1$t <- t
    data1$s <- s
    
    # fit FDboost function
    m2 <- FDboost(Y ~ 1 + bhist(X1, 
                                time = t, 
                                s = s, 
                                df = 5, 
                                knots = 25, 
                                limits = lD),
                  control = boost_control(mstop = 1000),
                  timeformula = ~ bbs(t),
                  data = data1)

    ppt <- cvrisk(m2, papply = mclapply, folds = cvLong(m2$id, type="kfold",
                                                        B = 15), mc.cores = coresCV)
    
    g1 <- coef(m2[mstop(ppt)],which=2)
    g1c <- g1$smterms[[1]]$value
    g1c[lower.tri(g1c,diag=T)] <- NA
    
    pred <- p[[2]]

    g2 <- matrix(pred$fit,ncol=obsPerTra)
    g2[lower.tri(g2,diag=T)] <- NA

    # calculate the reliMSEs of both approaches
    resDF[[nrSim]] <- data.frame(reliMSE_FDboost=sum(c(((g1c-trueBeta)^2)),na.rm = T) / 
                                   sum(c(trueBeta^2),na.rm = T),
                                 reliMSE_pffr=sum(c(((g2-trueBeta)^2)),na.rm = T) / 
                                   sum(c(trueBeta^2),na.rm = T)
    )
    
  }
  
  return(resDF)
  
}, mc.cores = coresSettings)

# save results
saveRDS(res,file="results/mm_sim2.RDS")

##########################################################################################
##########################################################################################

### do the same simulation with a different limit function

# new limit function -> used for band structure simulation
lD <- function(s, t) s < t & s >= t-0.10 & t <= 0.75

# parallelize simulations over different settings

res <- mclapply(1:nrow(setupDF),function(set){
  
  nn = setupDF[set,]$nn
  ss = setupDF[set,]$ss
  bd = setupDF[set,]$bd
  k = setupDF[set,]$k
  
  resDF <- vector("list",nrSims)
  
  for(nrSim in 1:nrSims){

    # 88, 77, 18
        
    obsPerTra = 40
    data1 <- pffrSimVar(scenario = c("int", "ff"), 
                        n = nn, 
                        nxgrid = obsPerTra, 
                        nygrid = obsPerTra, 
                        limits = lD, 
                        SNR = ss,
                        simFun = hillyFun, 
                        bd = bd)
    
    
    data1$X1 <- t(t(data1$X1) - colMeans(data1$X1))
    
    X1 <- matrix(data1$X1,ncol=nn)
    Y <- matrix(data1$Y,ncol=nn)
    
    test1 <- hillyFun
    t <- attr(data1, "yindex")
    s <- attr(data1, "xindex")
    
    trueBeta <- outer(s, t, test1) * outer(s, t, lD)
    
    trueBeta[trueBeta==0] <- NA
    
    m1 <- pffr(Y ~ 1 + ff(X1, 
                          xind = s, 
                          splinepars = list(bs = "ps", m = list(2,1), k = k), 
                          yind = t,
                          limits = lD), 
               data = data1)
    p <- plot(m1, pers = FALSE, select = 2, scheme = 2, 
              rug = FALSE, n2 = obsPerTra)
    
    data1 <- as.list(data1)
    data1$t <- t
    data1$s <- s
    
    pred <- p[[2]]
    
    g2 <- matrix(pred$fit,ncol=obsPerTra)
    g2[lower.tri(g2,diag=T)] <- NA
    
    reliMSE_FDboost <- NA
    
    m2 <- try(FDboost(Y ~ 1 + bhist(X1, 
                                    time = t, 
                                    s = s, 
                                    df = 5, 
                                    knots = 25, 
                                    limits = lD),
                  control = boost_control(mstop = 1000),
                  timeformula = ~ bbs(t),
                  data = data1))

    if(class(m2)!="try-error"){
      
      ppt <- try(cvrisk(m2, papply = mclapply, 
                        folds = cvLong(m2$id, type = "kfold",
                                       B = 15), mc.cores = coresCV))
      
      if(class(ppt)!="try-error"){
        
        g1 <- coef(m2[mstop(ppt)],which=2)
        g1c <- g1$smterms[[1]]$value
        g1c[lower.tri(g1c,diag=T)] <- NA
        g1c[abs(g1c)<1e-16] <- NA

        reliMSE_FDboost = sum(c(((g1c-trueBeta)^2)),na.rm = T)/sum(c(trueBeta^2),na.rm = T)
        
      }
    }
    
    resDF[[nrSim]] <- data.frame(reliMSE_FDboost=reliMSE_FDboost,
                                 reliMSE_pffr=sum(c(((g2-trueBeta)^2)),na.rm = T)/sum(c(trueBeta^2),na.rm = T)
    )
    
  }
  
  return(resDF)
  
}, mc.cores = coresSettings)

# save results
saveRDS(res,file="results/band_sim2.RDS")
