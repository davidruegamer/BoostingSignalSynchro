############################################################################
### Simulation code for section 4.3 (uncertainty quantification)
############################################################################

source("0_libs_funs.R", chdir = T)
if(length(list.files("results")) == 0) dir.create("results")

nrSims = 100
## if you just want to test the code:
if(FALSE) nrSims = 2

### core usage
coresBoot = 25


# fix setting
obsPerTra <- c(40)
nuC <- 0.1
setup = "histOnly"
nrSims = 100
bootNr = 100
n <- c(192)
SNR <- c(0.1)

set.seed(42)

dat <- dataGenProc(n = n,
                   obsPerTra = obsPerTra,
                   SNR = SNR,
                   seed = 12,
                   setup = setup,
                   nrOfResp = nrSims, 
                   trueEffHist = function(s,t){ # true historical effect
                     
                     ret <- outer(s, t, function(s, t) #sin(abs(s-t))*cos(10*t)) # 
                       sin(abs(t-s)+10)*cos(5*s))
                          #sin(3*abs(t-s)+2)*cos(2*s-0.75))
                            #sin(abs(t-s+.5)+2.5))
                            # dnorm(s,0.25,0.05)*dnorm(t,0.6,0.05)/30)
                            ret[ret<1e-3] <- 0
                     
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

################################ model fit

# do for nrSims iterations
for(nrSim in 1:nrSims){
  
  set.seed(nrSim)
  
  dat$Yi <- if(nrSims==1) dat$Y else  dat$Y[[nrSim]]
  dat$Yvec <- Yvec <- as.vector(dat$Yi)
  
  # fit initial model
  mod2 <- FDboost(Yi ~ 1 + bhistx(X1h, df = 15, knots = 5, differences = 1, standard = 'length'),
                  timeformula = ~ bbs(t, df=2.5), 
                  data = dat,
                  control = boost_control(mstop = 500, nu = nuC)
  )
  
  # bootstrap model
  bootR <- try(bootMod(mod = mod2, 
                       nrBoots = bootNr, 
                       which = 2,
                       data = dat, 
                       obsPerTra = obsPerTra, 
                       vars = c("g2","g3","X1h","Yi"),
                       argvals = c("s","t"), 
                       idvars = c("repIDx"),
                       mcCores = coresBoot))
 
  saveRDS(bootR,file=paste0("results/bootHistOnly4/boot_histOnly_nrSim",nrSim,".RDS"))
  
}
