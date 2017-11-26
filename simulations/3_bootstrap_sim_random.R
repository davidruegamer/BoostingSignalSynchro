############################################################################
### Simulation code for section 4.3 (uncertainty quantification)
############################################################################

source("0_libs_funs.R", chdir = T)
if(length(list.files("results")) == 0) dir.create("results")
if(length(list.files("results/bootRandom2")) == 0) dir.create("results/bootRandom2")


nrSims = 100
## if you just want to test the code:
if(FALSE) nrSims = 2

### core usage
coresBoot = 25

# fix setting
n <- c(160
       )
SNR <- c(1)
obsPerTra <- c(#40, 
  40)
nuC <- 0.1
setup = "histRandIA"
bootNr = 100
## if you just want to test the code:
if(FALSE) bootNr = 2
nrRanEf = 10

set.seed(42)

# generate data
dat <- dataGenProc(n = n,
                   obsPerTra = obsPerTra,
                   SNR = SNR,
                   seed = 12,
                   setup = setup,
                   nrRanEf = nrRanEf,
                   nrOfResp = nrSims, 
                   trueEffHist = function(s,t){ # true historical effect
                     
                     ret <- outer(s, t, function(s, t) sin(abs(s-t))*cos(10*t)*15)
                     ret[ret<1e-4] <- 0
                     
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

# do for 100 iterations
mclapply(1:nrSims, function(nrSim){
  
  set.seed(nrSim)
  
  dat$Yi <- if(nrSims == 1) dat$Y else  dat$Y[[nrSim]]
  dat$Yvec <- Yvec <- as.vector(dat$Yi)
  
  # fit initial model
  mod2 <- FDboost(Yi ~ 1 + brandomc(g3, df = 6) +
                    bhistx(X1h, df = 15, knots = 5, differences = 1, standard = 'length') + 
                  bhistx(X1h, df = 6, knots = 5, differences = 1, standard = 'length') %X% 
                    brandomc(g3, index = repIDx, df = 2.5),
                  timeformula = ~ bbs(t, df = 2.5), data = dat,
                  control = boost_control(mstop = 5000, nu = nuC)
  )
  
  
  by <- "g3"

  print("Done with fit. Bootstrapping...")

  # bootstrap models
  bootR <- try(bootMod(mod = mod2, 
                       nrBoots = bootNr, 
                       which=3, 
                       vars = c("X1h","Yi","g2","g3"), 
                       argvals = c("s","t"),
                       idvars = c("repIDx"), 
                       retBootCoefs=TRUE,
                       data = dat, 
                       obsPerTra = obsPerTra, 
                       by = by,
                       mcCores = coresBoot))
  
  saveRDS(bootR, file = paste0("results/bootRandom2/boot_randIA_nrSim", nrSim, ".RDS"))
  
}, mc.cores=1)
