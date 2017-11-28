library(coda)

source("simpleExpo.R")
source("simpleWeibull.R")
source("simpleBAOWei.R")
source("simplePLP.R")
source("simpleMPLP.R")

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
# Expo test
################################################################################

expoResults <- simpExpoMCMC(subListGAN[[20]], 10000, burnin = 2000)

expoResults$DIC
expoResults$PD

plot(as.mcmc(expoResults$lam_draws))
plot(as.mcmc(expoResults$theta1_draws))
plot(as.mcmc(expoResults$theta2_draws))


################################################################################
# Weibull test
################################################################################

system.time(weiResults <- simpWeiMCMC(subListGAN[[20]], 10000, burnin = 2000))

weiResults$DIC
weiResults$PD

plot(as.mcmc(weiResults$lam_draws))
plot(as.mcmc(weiResults$theta1_draws))
plot(as.mcmc(weiResults$theta2_draws))
plot(as.mcmc(weiResults$shape_draws))


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
# Power Law Process test
################################################################################
PLPResults1 <- simpPLP(subList[[20]], 5000, burnin = 2000)

PLPResults1$DIC
PLPResults1$PD

plot(as.mcmc(PLPResults1$draws[,3]))
plot(as.mcmc(PLPResults1$draws[,4]))
plot(as.mcmc(PLPResults1$draws[,5]))
plot(as.mcmc(PLPResults1$draws[,6]))

################################################################################
# Weibull test
################################################################################

weiResults2 <- simpBAOWeiMCMC(subList[[20]], 20000, burnin = 5000)

weiResults2$DIC
weiResults2$PD

plot(as.mcmc(weiResults2$lam_draws))
plot(as.mcmc(weiResults2$theta1_draws))
plot(as.mcmc(weiResults2$theta2_draws))
plot(as.mcmc(weiResults2$shape_draws))


###########################################################################

# simulation to test code MPLP
shape <- .5
lambda <- rep(1, 8)
theta1 <- .8
theta2 <- .9
kappa <- 1.5

L.fun <- function(t, scale, shape){
  return(scale * t^shape)
}
L.inv <- function(d, scale, shape){
  return((d/scale)^(1/shape))
}

testData3 <- list()

for(j in 1:26){
  testData3[[j]] <- vector("list", 8)
  
  for(i in 1:8){
    testData3[[j]][[i]] <- data.frame(totalMiles = rep(0, 1000),
                                      Censor = c(rep(0, 399), 1, rep(0, 299), 1, rep(0, 299), 1),
                                      Phase = c(rep(1,400), rep(2,300), rep(3,300)),
                                      trun = rep(F, 1000),
                                      Lower = rep(0, 1000),
                                      Upper = rep(0, 1000),
                                      MBF = rep(0, 1000))
    testData3[[j]][[i]]$totalMiles[1] <- L.inv(rgamma(1, kappa, 1), lambda[i], shape)
    for(k in 2:400){
      testData3[[j]][[i]]$totalMiles[k] <- L.inv(rgamma(1, kappa, 1) + L.fun(testData3[[j]][[i]]$totalMiles[k-1], lambda[i], shape), lambda[i], shape)
    }
    testData3[[j]][[i]]$totalMiles[401] <- L.inv(rgamma(1, kappa, 1), lambda[i]*theta1, shape)
    for(k in 402:700){
      testData3[[j]][[i]]$totalMiles[k] <- L.inv(rgamma(1, kappa, 1) + L.fun(testData3[[j]][[i]]$totalMiles[k-1], lambda[i]*theta1, shape), lambda[i]*theta1, shape)
    }
    testData3[[j]][[i]]$totalMiles[701] <- L.inv(rgamma(1, kappa, 1), lambda[i]*theta1*theta2, shape)
    for(k in 702:1000){
      testData3[[j]][[i]]$totalMiles[k] <- L.inv(rgamma(1, kappa, 1) + L.fun(testData3[[j]][[i]]$totalMiles[k-1], lambda[i]*theta1*theta2, shape), lambda[i]*theta1*theta2, shape)
    }
    testData3[[j]][[i]]$MBF[1] <- testData3[[j]][[i]]$totalMiles[1]
    testData3[[j]][[i]]$MBF[401] <- testData3[[j]][[i]]$totalMiles[401]
    testData3[[j]][[i]]$MBF[701] <- testData3[[j]][[i]]$totalMiles[701]
    testData3[[j]][[i]]$Upper[1] <- testData3[[j]][[i]]$totalMiles[1]
    testData3[[j]][[i]]$Upper[401] <- testData3[[j]][[i]]$totalMiles[401]
    testData3[[j]][[i]]$Upper[701] <- testData3[[j]][[i]]$totalMiles[701]
    for(k in 2:400){
      testData3[[j]][[i]]$MBF[k] <- testData3[[j]][[i]]$totalMiles[k] - testData3[[j]][[i]]$totalMiles[k-1]
      testData3[[j]][[i]]$Lower[k] <- testData3[[j]][[i]]$totalMiles[k-1]
      testData3[[j]][[i]]$Upper[k] <- testData3[[j]][[i]]$totalMiles[k]
    }
    for(k in 402:700){
      testData3[[j]][[i]]$MBF[k] <- testData3[[j]][[i]]$totalMiles[k] - testData3[[j]][[i]]$totalMiles[k-1]
      testData3[[j]][[i]]$Lower[k] <- testData3[[j]][[i]]$totalMiles[k-1]
      testData3[[j]][[i]]$Upper[k] <- testData3[[j]][[i]]$totalMiles[k]
    }
    for(k in 702:1000){
      testData3[[j]][[i]]$MBF[k] <- testData3[[j]][[i]]$totalMiles[k] - testData3[[j]][[i]]$totalMiles[k-1]
      testData3[[j]][[i]]$Lower[k] <- testData3[[j]][[i]]$totalMiles[k-1]
      testData3[[j]][[i]]$Upper[k] <- testData3[[j]][[i]]$totalMiles[k]
    }
  }
}

source("simpleMPLP.R")
MPLPResults <- simpMPLPMCMC(testData3[[1]], 20000, burnin = 5000)

MPLPResults$DIC
MPLPResults$PD

print(lambda)

plot(as.mcmc(MPLPResults$draws[,3]))
plot(as.mcmc(MPLPResults$draws[,4]))
plot(as.mcmc(MPLPResults$draws[,5]))
plot(as.mcmc(MPLPResults$draws[,6]))
plot(as.mcmc(MPLPResults$draws[,7]))



