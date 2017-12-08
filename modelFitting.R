
library(compiler)
library(parallel)

################################################################################
# JLTV Data Good As New testing
################################################################################
JLTVData <- read.csv("JLTVData.csv")

JLTVData[, 1] <- as.numeric(JLTVData[, 1])
JLTVData$trun <- JLTVData$MBF == 0

JLTVData$Lower <- 0
JLTVData$Upper <- 0

subs <- list()

for(i in 1:26){
  subs[[i]] <- JLTVData[which(JLTVData[,9] == i),]
}

subListGAN <- list()

for(i in 1:26){
  subListGAN[[i]] <- list()
  for(j in 1:8){
    subListGAN[[i]][[j]] <- subs[[i]][which(subs[[i]][,1]==j),c(6,8,3,10:12)]
  }
}

for(i in 1:length(subListGAN)){
  for(j in 1:length(subListGAN[[i]]))
    if(all(subListGAN[[i]][[j]]$trun == F)){
    } else {
      for(k in 1:nrow(subListGAN[[i]][[j]]))
        if(subListGAN[[i]][[j]][k,4] == T){
          if(k == 1){
          } else {
            
            if(subListGAN[[i]][[j]][k-1,2] == 1 & subListGAN[[i]][[j]][k-1,3] != 3){
              subListGAN[[i]][[j]][k,3] <- subListGAN[[i]][[j]][k-1,3] + 1
            }
            
            if(subListGAN[[i]][[j]][k+1,4] == F){
              subListGAN[[i]][[j]][k,6] <- subListGAN[[i]][[j]][k+1,1]
            } else if(subListGAN[[i]][[j]][k+2,4] == F){
              subListGAN[[i]][[j]][k,6] <- subListGAN[[i]][[j]][k+2,1]
            } else if(subListGAN[[i]][[j]][k+3,4] == F){
              subListGAN[[i]][[j]][k,6] <- subListGAN[[i]][[j]][k+3,1]
            } else if(subListGAN[[i]][[j]][k+4,4] == F){
              subListGAN[[i]][[j]][k,6] <- subListGAN[[i]][[j]][k+4,1]
            } else if(subListGAN[[i]][[j]][k+5,4] == F){
              subListGAN[[i]][[j]][k,6] <- subListGAN[[i]][[j]][k+5,1]
            } else{
              subListGAN[[i]][[j]][k,6] <- subListGAN[[i]][[j]][k+6,1]
            }
          }
        }
    }
}

################################################################################
# JLTV Data Bad As Old testing
################################################################################

JLTVData <- read.csv("JLTVData.csv")

JLTVData[, 1] <- as.numeric(JLTVData[, 1])
JLTVData$trun <- JLTVData$MBF == 0

JLTVData$Lower <- 0
JLTVData$Upper <- 0

JLTVData$totalMiles <- 0

subs <- list()

for(i in 1:26){
  subs[[i]] <- JLTVData[which(JLTVData[,9] == i),]
}

subList <- list()

for(i in 1:26){
  subList[[i]] <- list()
  for(j in 1:8){
    subList[[i]][[j]] <- subs[[i]][which(subs[[i]][,1]==j),c(6,8,3,10:13)]
  }
}


for(i in 1:length(subList)){
  for(j in 1:length(subList[[i]])){
    for(k in 1:3){
      failIndex <- which(subList[[i]][[j]][,3] == k)
      subList[[i]][[j]][failIndex[1],7] <- subList[[i]][[j]][failIndex[1],1]
      if(length(failIndex) > 1){
        for(l in failIndex[-1]){
          subList[[i]][[j]][l,7] <- subList[[i]][[j]][l,1] + subList[[i]][[j]][l-1,7]
        }
      }
    }
  }     
}

for(i in 1:length(subList)){
  for(j in 1:length(subList[[i]]))
    for(k in 1:nrow(subList[[i]][[j]]))
      if(subList[[i]][[j]][k,4] == T){
        if(k == 1 || subList[[i]][[j]][k,3] > subList[[i]][[j]][k-1,3]){
          subList[[i]][[j]][k,6] <- subList[[i]][[j]][k,7]
        } else {
          if(subList[[i]][[j]][k-1,2] == 1 & subList[[i]][[j]][k-1,3] != 3){
            subList[[i]][[j]][k,3] <- subList[[i]][[j]][k-1,3] + 1
          }
          if(subList[[i]][[j]][k,3] != subList[[i]][[j]][k-1,3]){
            subList[[i]][[j]][k,7] <- 0
          }
          subList[[i]][[j]][k,5] <- subList[[i]][[j]][k,7]
          if(subList[[i]][[j]][k+1,4] == F){
            subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+1,7]
          } else if(subList[[i]][[j]][k+2,4] == F){
            subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+2,7]
          } else if(subList[[i]][[j]][k+3,4] == F){
            subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+3,7]
          } else if(subList[[i]][[j]][k+4,4] == F){
            subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+4,7]
          } else if(subList[[i]][[j]][k+5,4] == F){
            subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+5,7]
          } else{
            subList[[i]][[j]][k,6] <- subList[[i]][[j]][k+6,7]
          }
        }
      } else {
        if(k == 1 || subList[[i]][[j]][k,3] > subList[[i]][[j]][k-1,3]){
          subList[[i]][[j]][k,6] <- subList[[i]][[j]][k,7]
        } else {
          subList[[i]][[j]][k,6] <- subList[[i]][[j]][k,7]
          subList[[i]][[j]][k,5] <- subList[[i]][[j]][k-1,7]
        }
      }
}
################################################################################
# Running simple models
################################################################################

ncores <- detectCores() - 2

source("simpleExpo.R")
M1 <- cmpfun(simpExpoMCMC)

cl <- makeCluster(ncores)
system.time(simpExpoResults <- parLapply(cl,subListGAN, M1))
stopCluster(cl)
save(simpExpoResults ,file="simpExpoResults.RData")

source("simpleWeibull.R")
M2 <- cmpfun(simpWeiMCMC)

cl <- makeCluster(ncores)
clusterExport(cl,c("logWeibullSum"))
system.time(simpWeiResults <- parLapply(cl,subListGAN, M2))
stopCluster(cl)
save(simpWeiResults ,file="simpWeiResults.RData")

source("simpleBAOWei.R")
M3 <- cmpfun(simpBAOWeiMCMC)

cl <- makeCluster(ncores)
clusterExport(cl,c("logWeibullSum"))
system.time(simpBAOWeiResults <- parLapply(cl,subList, M3))
stopCluster(cl)
save(simpBAOWeiResults ,file="simpBAOWeiResults.RData")

source("simplePLP.R")
M4 <- cmpfun(simpPLP)

cl <- makeCluster(ncores)
clusterExport(cl,c("logPLPSum","logPLPInterval"))
system.time(simpPLPResults <- parLapply(cl,subList, M4))
stopCluster(cl)
save(simpPLPResults ,file="simpPLPResults.RData")

source("simpleMPLP.R")
M5 <- cmpfun(simpMPLPMCMC)

cl <- makeCluster(ncores)
clusterExport(cl,c("logMPLPsum", "logMPLPInterval"))
system.time(simpMPLPResults <- parLapply(cl,subList, M5))
stopCluster(cl)
save(simpMPLPResults ,file="simpMPLPResults.RData")

################################################################################
# Running hierarchical models
################################################################################

source("expoMCMC.R")
M6 <- cmpfun(expoMCMC)

cl <- makeCluster(ncores)
system.time(ExpoResults <- parLapply(cl,subListGAN, M6))
stopCluster(cl)
save(ExpoResults ,file="ExpoResults.RData")

source("weibullMCMC.R")
M7 <- cmpfun(WeibullMCMC)

cl <- makeCluster(ncores)
clusterExport(cl,c("logWeibullSum"))
system.time(WeibullResults <- parLapply(cl,subListGAN, M7))
stopCluster(cl)
save(WeibullResults ,file="WeibullResults.RData")

source("BAOWeibullMCMC.R")
M8 <- cmpfun(BAOWMCMC)

cl <- makeCluster(ncores)
clusterExport(cl,c("logWeibullSum"))
system.time(BAOWeiResults <- parLapply(cl,subList, M8))
stopCluster(cl)
save(BAOWeiResults ,file="BAOWeiResults.RData")

source("PLPMCMC.R")
M9 <- cmpfun(PLPMCMC)

cl <- makeCluster(ncores)
clusterExport(cl,c("logPLPSum","logPLPInterval"))
system.time(PLPResults <- parLapply(cl,subList, M9))
stopCluster(cl)
save(PLPResults ,file="PLPResults.RData")

source("MPLPMCMC.R")
M10 <- cmpfun(MPLPMCMC)

cl <- makeCluster(ncores)
clusterExport(cl,c("logMPLPsum", "logMPLPInterval"))
system.time(MPLPResults <- parLapply(cl,subList, M10))
stopCluster(cl)
save(MPLPResults ,file="MPLPResults.RData")


DICResults <- matrix(0, nrow = 26, ncol = 10)
for(i in 1:26){
  DICResults[i,1] <- simpExpoResults[[i]]$DIC
  DICResults[i,2] <- simpWeiResults[[i]]$DIC
  DICResults[i,3] <- simpBAOWeiResults[[i]]$DIC
  DICResults[i,4] <- simpPLPResults[[i]]$DIC
  DICResults[i,5] <- simpMPLPResults[[i]]$DIC

  DICResults[i,6] <- ExpoResults[[i]]$DIC
  DICResults[i,7] <- WeibullResults[[i]]$DIC
  DICResults[i,8] <- BAOWeiResults[[i]]$DIC
  DICResults[i,9] <- PLPResults[[i]]$DIC
  DICResults[i,10] <- MPLPResults[[i]]$DIC
}

print(DICResults)

tot <- 0
for(i in 1:26){
  if(!is.na(DICResults[i,5])){
    if(is.finite(DICResults[i,5])){
      tot <- tot + DICResults[i,5]
    } else {
      tot <-   tot + DICResults[i,4]
    }
  } else {
    tot <- tot + DICResults[i,2]
  }
}
tot

tot <- 0
for(i in 1:21){
  tot <- tot + min(DICResults[i, 3:5])
}
tot <- tot + sum(DICResults[c(13, 15, 17, 18, 23,24, 26), 4]) + sum(DICResults[c(14, 16, 19:21, 25), 5]) + DICResults[22, 2]

tot

tot <- 0
for(i in c(1:21, 24:26)){
  tot <- tot + min(DICResults[i, c(9,10)])
}
tot <- tot + DICResults[22,7] + DICResults[23,9]
tot

tot <- 0
for(i in c(1:12, 14:16, 19:21, 24:25)){
  tot <- tot + min(DICResults[i, c(3:5, 9,10)])
}
tot <- tot + DICResults[13,10] + DICResults[17,9] + DICResults[18,4] + DICResults[22,7] + DICResults[23,9] + DICResults[26,9]
tot

sum(DICResults[-7,8])
sum(DICResults[c(-10, -22),3])

tot <- 0
for(i in 1:26){
  tot <- tot + min(DICResults[i,c(3,8)])
}
tot

for(i in 1:26){
  print(DICResults[i,8] - DICResults[i,7])
}

for(i in 1:26){
  print(DICResults[i,3] - DICResults[i,2])
}

for(i in 1:26){
  print(DICResults[i,8] - DICResults[i,3])
}

sum(DICResults[1:2,8], DICResults[3,3], DICResults[4:6,8], DICResults[7,3], DICResults[8:26,8])

plot(as.mcmc(MPLPResults[[20]]$draws[,1]))
autocorr.plot(MPLPResults[[20]]$draws[,1])

plot(as.mcmc(MPLPResults[[20]]$draws[,2]))
autocorr.plot(MPLPResults[[20]]$draws[,2])

plot(as.mcmc(MPLPResults[[20]]$draws[,3]))
autocorr.plot(MPLPResults[[20]]$draws[,3])

plot(as.mcmc(MPLPResults[[20]]$draws[,4]))
autocorr.plot(MPLPResults[[20]]$draws[,4])

plot(as.mcmc(MPLPResults[[20]]$draws[,5]))
autocorr.plot(MPLPResults[[20]]$draws[,5])

plot(as.mcmc(MPLPResults[[20]]$draws[,6]))
autocorr.plot(MPLPResults[[20]]$draws[,6])

plot(as.mcmc(MPLPResults[[20]]$draws[,7]))
autocorr.plot(MPLPResults[[20]]$draws[,7])

plot(as.mcmc(MPLPResults[[20]]$draws[,8]))
autocorr.plot(MPLPResults[[20]]$draws[,8])

plot(as.mcmc(MPLPResults[[20]]$draws[,9]))
autocorr.plot(MPLPResults[[20]]$draws[,9])

plot(as.mcmc(MPLPResults[[20]]$draws[,10]))
autocorr.plot(MPLPResults[[20]]$draws[,10])

plot(as.mcmc(MPLPResults[[20]]$draws[,11]))
autocorr.plot(MPLPResults[[20]]$draws[,11])

plot(as.mcmc(MPLPResults[[20]]$draws[,12]))
autocorr.plot(MPLPResults[[20]]$draws[,12])

plot(as.mcmc(MPLPResults[[20]]$draws[,13]))
autocorr.plot(MPLPResults[[20]]$draws[,13])

plot(as.mcmc(MPLPResults[[20]]$draws[,14]))
autocorr.plot(MPLPResults[[20]]$draws[,14])


plot(as.mcmc(simpExpoResults[[3]]$lam_draws))
plot(as.mcmc(simpExpoResults[[3]]$theta1_draws))
plot(as.mcmc(simpExpoResults[[3]]$theta2_draws))


plot(as.mcmc(simpBAOWeiResults[[10]]$lam_draws))
plot(as.mcmc(simpBAOWeiResults[[10]]$theta1_draws))
plot(as.mcmc(simpBAOWeiResults[[10]]$theta2_draws))
plot(as.mcmc(simpBAOWeiResults[[10]]$shape_draws))

plot(as.mcmc(simpBAOWeiResults[[3]]$lam_draws))
plot(as.mcmc(simpBAOWeiResults[[3]]$theta1_draws))
plot(as.mcmc(simpBAOWeiResults[[3]]$theta2_draws))
plot(as.mcmc(simpBAOWeiResults[[3]]$shape_draws))

simpBAOWeiResults[[3]]$DIC
simpBAOWeiResults[[3]]$PD

plot(as.mcmc(WeibullResults[[3]]$draws[,1]))
plot(as.mcmc(WeibullResults[[3]]$draws[,2]))
plot(as.mcmc(WeibullResults[[3]]$draws[,3]))
plot(as.mcmc(WeibullResults[[3]]$draws[,4]))
plot(as.mcmc(WeibullResults[[3]]$draws[,5]))
plot(as.mcmc(WeibullResults[[3]]$draws[,6]))

plot(as.mcmc(ExpoResults[[3]]$draws[,1]))
plot(as.mcmc(ExpoResults[[3]]$draws[,2]))
plot(as.mcmc(ExpoResults[[3]]$draws[,3]))
plot(as.mcmc(ExpoResults[[3]]$draws[,4]))
plot(as.mcmc(ExpoResults[[3]]$draws[,5]))
plot(as.mcmc(ExpoResults[[3]]$draws[,6]))

ExpoResults[[3]]$acceptance
ExpoResults[[3]]$PD

plot(as.mcmc(simpMPLPResults[[1]]$draws[,3]))
plot(as.mcmc(simpMPLPResults[[1]]$draws[,4]))
plot(as.mcmc(simpMPLPResults[[1]]$draws[,5]))
plot(as.mcmc(simpMPLPResults[[1]]$draws[,6]))
plot(as.mcmc(simpMPLPResults[[1]]$draws[,7]))

plot(as.mcmc(simpMPLPResults[[4]]$draws[,3]))
plot(as.mcmc(simpMPLPResults[[4]]$draws[,4]))
plot(as.mcmc(simpMPLPResults[[4]]$draws[,5]))
plot(as.mcmc(simpMPLPResults[[4]]$draws[,6]))
plot(as.mcmc(simpMPLPResults[[4]]$draws[,7]))

simpMPLPResults[[24]]$acceptanceKappa
simpMPLPResults[[24]]$acceptanceShape


plot(as.mcmc(BAOWeiResults[[7]]$draws[,4]))
plot(as.mcmc(BAOWeiResults[[7]]$draws[,5]))
plot(as.mcmc(BAOWeiResults[[7]]$draws[,6]))
plot(as.mcmc(BAOWeiResults[[7]]$draws[,7]))
plot(as.mcmc(BAOWeiResults[[7]]$draws[,8]))
plot(as.mcmc(BAOWeiResults[[7]]$draws[,9]))

for(i in 1:26){
  print(which(DICResults[i,] == min(DICResults[i,c(1:4, 6:9)])))
}


plot(density(WeibullResults[[1]]$draws[,3]), ylim = c(0,50), xlim = c(0,1.5), 
     main = "Posterior Plots of Shape Parameter", xlab = "")
for(i in 2:26){
  lines(density(WeibullResults[[i]]$draws[,3]))
}

plot(density(MPLPResults[[1]]$draws[,3]), ylim = c(0,50), xlim = c(0,5), 
     main = "Posterior Plots of Shock Parameter", xlab = "")
for(i in 2:26){
  lines(density(MPLPResults[[i]]$draws[,3]))
}