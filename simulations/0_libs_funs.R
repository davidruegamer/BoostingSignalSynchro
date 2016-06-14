### libraries and sourcing

# install.packages("FDboost", repos="http://R-Forge.R-project.org")
library("parallel")
library("FDboost")
library("splines")
library("fda")

source("../functions/ranHist.R")
source("../functions/createRandomRespFolds.R")
source("../functions/dataGenProc.R")
source("../functions/bootMod.R")
