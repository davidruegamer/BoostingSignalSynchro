############################################################################
### Simulation code for section 4.3 (uncertainty quantification)
############################################################################

source("0_libs_funs.R", chdir = T)
if(length(list.files("results")) == 0) dir.create("results")
if(length(list.files("results/bootFac")) == 0) dir.create("results/bootFac")

nrSims = 100
## if you just want to test the code:
if(FALSE) nrSims = 2

### core usage
coresBoot = 25

# fix setting
obsPerTra <- c(40)
nuC <- 0.1
setup = "histGameIA"
bootNr = 100
## if you just want to test the code:
if(FALSE) bootNr = 2
n <- c(160)
SNR <- c(1)

set.seed(42)

dat <- dataGenProc(n = n,
                   nrFacEf = 4,
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
                   },
                   trueEffHistGame = function(s,t){ # true historical effect

                     ret <- outer(s, t, function(s, t) dnorm(s,0.9,0.2)*dnorm(t,0.9,0.2))
                     ret[ret<0.001] <- 0

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

# do for nrSims iterations
for(nrSim in 1:nrSims){
  
  set.seed(nrSim)
  
  dat$Yi <- if(nrSims==1) dat$Y else  dat$Y[[nrSim]]
  dat$Yvec <- Yvec <- as.vector(dat$Yi)
  
  # fit initial model
  mod2 <- FDboost(Yi ~ 1 + bolsc(g2, df=1) +
                    bhistx(X1h, df=10, knots=5, differences=2, standard='length') + 
                  bhistx(X1h, df=10, knots=5, differences=2, standard='length') %X% 
                    bolsc(g2, index=repIDx, df=1),
                  timeformula = ~ bbs(t, df=10), data=dat,
                  control=boost_control(mstop=1500,nu=nuC)
  )
  
  # bootstrap model
  bootR <- try(bootMod(mod = mod2, 
                       nrBoots = bootNr, 
                       which = 3:4,
                       data = dat, 
                       obsPerTra = obsPerTra, 
                       vars = c("X1h","Yi","g2","g3"), 
                       argvals = c("s","t"),
                       idvars = c("repIDx"), 
                       retBootCoefs = TRUE,
                       mcCores = coresBoot))
  
  saveRDS(bootR,file=paste0("results/bootFac/boot_facIA_nrSim",nrSim,".RDS"))
  
}
