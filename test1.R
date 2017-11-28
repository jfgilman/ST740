library(coda)

################################################################################
# Testing Exponential-Gamma Hierarchical Model without censoring
################################################################################
source("expoMCMC.R")
source("weibullMCMC.R")
source("BAOWeibullMCMC.R")
source("PLPMCMC.R")
source("MPLPMCMC.R")
source("WeibullCDFfuns.R")

# simulation to test code
shape <- 1.5
lambda <- rep(c(.2,.5,1,3), 10)
theta1 <- .8
theta2 <- .9

testData <- list()

for(j in 1:26){
  testData[[j]] <- vector("list", 8)
  
  for(i in 1:8){
    u <- runif(1000)
    testData[[j]][[i]] <- data.frame(MBF = c(wF.inv(u[1:400],lambda[i],shape),
                                             wF.inv(u[401:700],lambda[i]*theta1,shape),
                                             wF.inv(u[701:1000],lambda[i]*theta1*theta2,shape)),
                                     Censor = rep(0, 1000),
                                     Phase = c(rep(1,400), rep(2,300), rep(3,300)),
                                     trun = rep(F, 1000))
  }
}

expoResults1 <- expoMCMC(testData[[1]], 10000, tuning = 8)

expoResults1$acceptance
expoResults1$DIC
expoResults1$PD

plot(as.mcmc(expoResults1$draws[,1]))
plot(as.mcmc(expoResults1$draws[,2]))
plot(as.mcmc(expoResults1$draws[,3]))
plot(as.mcmc(expoResults1$draws[,4]))
plot(as.mcmc(expoResults1$draws[,5]))
plot(as.mcmc(expoResults1$draws[,6]))
plot(as.mcmc(expoResults1$draws[,7]))
plot(as.mcmc(expoResults1$draws[,8]))
plot(as.mcmc(expoResults1$draws[,9]))
plot(as.mcmc(expoResults1$draws[,10]))

################################################################################
# Weibull testing
################################################################################
weibullResults1 <- WeibullMCMC(testData[[1]], 8000)

weibullResults1$DIC
weibullResults1$PD

plot(as.mcmc(weibullResults1$draws[,1]))
plot(as.mcmc(weibullResults1$draws[,2]))
plot(as.mcmc(weibullResults1$draws[,3]))
plot(as.mcmc(weibullResults1$draws[,4]))
plot(as.mcmc(weibullResults1$draws[,5]))
plot(as.mcmc(weibullResults1$draws[,6]))
plot(as.mcmc(weibullResults1$draws[,7]))
plot(as.mcmc(weibullResults1$draws[,8]))
plot(as.mcmc(weibullResults1$draws[,9]))
plot(as.mcmc(weibullResults1$draws[,10]))


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

expoResults2 <- expoMCMC(subListGAN[[20]], 10000, tuning = 8)

expoResults2$DIC
expoResults2$PD

plot(as.mcmc(expoResults2$draws[,1]))
plot(as.mcmc(expoResults2$draws[,2]))
plot(as.mcmc(expoResults2$draws[,3]))
plot(as.mcmc(expoResults2$draws[,4]))
plot(as.mcmc(expoResults2$draws[,5]))
plot(as.mcmc(expoResults2$draws[,6]))
plot(as.mcmc(expoResults2$draws[,7]))
plot(as.mcmc(expoResults2$draws[,8]))
plot(as.mcmc(expoResults2$draws[,9]))
plot(as.mcmc(expoResults2$draws[,10]))

################################################################################
# Weibull test
################################################################################

weibullResults2 <- WeibullMCMC(subListGAN[[20]], 10000)

weibullResults2$DIC
weibullResults2$PD

plot(as.mcmc(weibullResults2$draws[,1]))
plot(as.mcmc(weibullResults2$draws[,2]))
plot(as.mcmc(weibullResults2$draws[,3]))
plot(as.mcmc(weibullResults2$draws[,4]))
plot(as.mcmc(weibullResults2$draws[,5]))
plot(as.mcmc(weibullResults2$draws[,6]))
plot(as.mcmc(weibullResults2$draws[,7]))
plot(as.mcmc(weibullResults2$draws[,8]))
plot(as.mcmc(weibullResults2$draws[,9]))
plot(as.mcmc(weibullResults2$draws[,10]))



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
# Weibull test
################################################################################

weibullResults3 <- BAOWMCMC(subList[[20]], 20000, burnin = 5000)

weibullResults3$DIC
weibullResults3$PD

plot(as.mcmc(weibullResults3$draws[,1]))
plot(as.mcmc(weibullResults3$draws[,2]))
plot(as.mcmc(weibullResults3$draws[,3]))
plot(as.mcmc(weibullResults3$draws[,4]))
plot(as.mcmc(weibullResults3$draws[,5]))
plot(as.mcmc(weibullResults3$draws[,6]))
plot(as.mcmc(weibullResults3$draws[,7]))
plot(as.mcmc(weibullResults3$draws[,8]))
plot(as.mcmc(weibullResults3$draws[,9]))
plot(as.mcmc(weibullResults3$draws[,10]))


################################################################################
# Power Law Process test
################################################################################

PLPResults1 <- PLPMCMC(subList[[20]], 20000, burnin = 5000)

PLPResults1$DIC
PLPResults1$PD

plot(as.mcmc(PLPResults1$draws[,1]))
plot(as.mcmc(PLPResults1$draws[,2]))
plot(as.mcmc(PLPResults1$draws[,3]))
plot(as.mcmc(PLPResults1$draws[,4]))
plot(as.mcmc(PLPResults1$draws[,5]))
plot(as.mcmc(PLPResults1$draws[,6]))
plot(as.mcmc(PLPResults1$draws[,7]))
plot(as.mcmc(PLPResults1$draws[,8]))
plot(as.mcmc(PLPResults1$draws[,9]))
plot(as.mcmc(PLPResults1$draws[,10]))

################################################################################
# Modulated Power Law Process test
################################################################################

MPLPResults1 <- MPLPMCMC(subList[[20]], 20000, burnin = 5000)

MPLPResults1$DIC
MPLPResults1$PD

plot(as.mcmc(MPLPResults1$draws[,1]))
plot(as.mcmc(MPLPResults1$draws[,2]))
plot(as.mcmc(MPLPResults1$draws[,3]))
plot(as.mcmc(MPLPResults1$draws[,4]))
plot(as.mcmc(MPLPResults1$draws[,5]))
plot(as.mcmc(MPLPResults1$draws[,6]))
plot(as.mcmc(MPLPResults1$draws[,7]))
plot(as.mcmc(MPLPResults1$draws[,8]))
plot(as.mcmc(MPLPResults1$draws[,9]))
plot(as.mcmc(MPLPResults1$draws[,10]))

# simulation to test code PLP
shape <- .3
lambda <- rgamma(8, .5, 1.5)
theta1 <- .8
theta2 <- .9

testData2 <- list()

for(j in 1:26){
  testData2[[j]] <- vector("list", 8)
  
  for(i in 1:8){
    testData2[[j]][[i]] <- data.frame(totalMiles = rep(0, 1000),
                                     Censor = c(rep(0, 399), 1, rep(0, 299), 1, rep(0, 299), 1),
                                     Phase = c(rep(1,400), rep(2,300), rep(3,300)),
                                     trun = rep(F, 1000),
                                     Lower = rep(0, 1000),
                                     Upper = rep(0, 1000),
                                     MBF = rep(0, 1000))
    testData2[[j]][[i]]$totalMiles[1] <- plpF.inv(runif(1), 0, lambda[i], shape)
    for(k in 2:400){
      testData2[[j]][[i]]$totalMiles[k] <- plpF.inv(runif(1), testData2[[j]][[i]]$totalMiles[k-1], lambda[i], shape)
    }
    testData2[[j]][[i]]$totalMiles[401] <- plpF.inv(runif(1), 0, lambda[i]*theta1, shape)
    for(k in 402:700){
      testData2[[j]][[i]]$totalMiles[k] <- plpF.inv(runif(1), testData2[[j]][[i]]$totalMiles[k-1], lambda[i]*theta1, shape)
    }
    testData2[[j]][[i]]$totalMiles[701] <- plpF.inv(runif(1), 0, lambda[i]*theta1*theta2, shape)
    for(k in 702:1000){
      testData2[[j]][[i]]$totalMiles[k] <- plpF.inv(runif(1), testData2[[j]][[i]]$totalMiles[k-1], lambda[i]*theta1*theta2, shape)
    }
    testData2[[j]][[i]]$MBF[1] <- testData2[[j]][[i]]$totalMiles[1]
    testData2[[j]][[i]]$MBF[401] <- testData2[[j]][[i]]$totalMiles[401]
    testData2[[j]][[i]]$MBF[701] <- testData2[[j]][[i]]$totalMiles[701]
    testData2[[j]][[i]]$Upper[1] <- testData2[[j]][[i]]$totalMiles[1]
    testData2[[j]][[i]]$Upper[401] <- testData2[[j]][[i]]$totalMiles[401]
    testData2[[j]][[i]]$Upper[701] <- testData2[[j]][[i]]$totalMiles[701]
    for(k in 2:400){
      testData2[[j]][[i]]$MBF[k] <- testData2[[j]][[i]]$totalMiles[k] - testData2[[j]][[i]]$totalMiles[k-1]
      testData2[[j]][[i]]$Lower[k] <- testData2[[j]][[i]]$totalMiles[k-1]
      testData2[[j]][[i]]$Upper[k] <- testData2[[j]][[i]]$totalMiles[k]
    }
    for(k in 402:700){
      testData2[[j]][[i]]$MBF[k] <- testData2[[j]][[i]]$totalMiles[k] - testData2[[j]][[i]]$totalMiles[k-1]
      testData2[[j]][[i]]$Lower[k] <- testData2[[j]][[i]]$totalMiles[k-1]
      testData2[[j]][[i]]$Upper[k] <- testData2[[j]][[i]]$totalMiles[k]
    }
    for(k in 702:1000){
      testData2[[j]][[i]]$MBF[k] <- testData2[[j]][[i]]$totalMiles[k] - testData2[[j]][[i]]$totalMiles[k-1]
      testData2[[j]][[i]]$Lower[k] <- testData2[[j]][[i]]$totalMiles[k-1]
      testData2[[j]][[i]]$Upper[k] <- testData2[[j]][[i]]$totalMiles[k]
    }
  }
}

PLPResults2 <- PLPMCMC(testData2[[1]], 10000, burnin = 5000)

PLPResults2$DIC
PLPResults2$PD

plot(as.mcmc(PLPResults2$draws[,1]))
plot(as.mcmc(PLPResults2$draws[,2]))
plot(as.mcmc(PLPResults2$draws[,3]))
plot(as.mcmc(PLPResults2$draws[,4]))
plot(as.mcmc(PLPResults2$draws[,5]))
plot(as.mcmc(PLPResults2$draws[,6]))
plot(as.mcmc(PLPResults2$draws[,7]))
plot(as.mcmc(PLPResults2$draws[,8]))
plot(as.mcmc(PLPResults2$draws[,9]))
plot(as.mcmc(PLPResults2$draws[,10]))



weibullResults4 <- BAOWMCMC(testData2[[1]], 10000, burnin = 5000)

weibullResults4$DIC
weibullResults4$PD

plot(as.mcmc(weibullResults4$draws[,1]))
plot(as.mcmc(weibullResults4$draws[,2]))
plot(as.mcmc(weibullResults4$draws[,3]))
plot(as.mcmc(weibullResults4$draws[,4]))
plot(as.mcmc(weibullResults4$draws[,5]))
plot(as.mcmc(weibullResults4$draws[,6]))
plot(as.mcmc(weibullResults4$draws[,7]))
plot(as.mcmc(weibullResults4$draws[,8]))
plot(as.mcmc(weibullResults4$draws[,9]))
plot(as.mcmc(weibullResults4$draws[,10]))


weibullResults5 <- WeibullMCMC(testData2[[1]], 10000, burnin = 5000)

weibullResults5$DIC
weibullResults5$PD

plot(as.mcmc(weibullResults5$draws[,1]))
plot(as.mcmc(weibullResults5$draws[,2]))
plot(as.mcmc(weibullResults5$draws[,3]))
plot(as.mcmc(weibullResults5$draws[,4]))
plot(as.mcmc(weibullResults5$draws[,5]))
plot(as.mcmc(weibullResults5$draws[,6]))
plot(as.mcmc(weibullResults5$draws[,7]))
plot(as.mcmc(weibullResults5$draws[,8]))
plot(as.mcmc(weibullResults5$draws[,9]))
plot(as.mcmc(weibullResults5$draws[,10]))



###########################################################################

# simulation to test code MPLP
shape <- .5
lambda <- rgamma(8, .5, 1.5)
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

MPLPResults2 <- MPLPMCMC(testData3[[1]], 20000, burnin = 10000)

MPLPResults2$DIC
MPLPResults2$PD

print(lambda)

plot(as.mcmc(MPLPResults2$draws[,1]))
plot(as.mcmc(MPLPResults2$draws[,2]))
plot(as.mcmc(MPLPResults2$draws[,3]))
plot(as.mcmc(MPLPResults2$draws[,4]))
plot(as.mcmc(MPLPResults2$draws[,5]))
plot(as.mcmc(MPLPResults2$draws[,6]))
plot(as.mcmc(MPLPResults2$draws[,7]))
plot(as.mcmc(MPLPResults2$draws[,8]))
plot(as.mcmc(MPLPResults2$draws[,9]))
plot(as.mcmc(MPLPResults2$draws[,10]))

