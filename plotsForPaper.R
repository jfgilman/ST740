

# Re run and do for rare data
load("BAOWeiResults.RData")
load("ExpoResults.RData")
load("MPLPResults.RData")
load("PLPResults.RData")
load("simpBAOWeiResults.RData")
load("simpExpoResults.RData")
load("simpMPLPResults.RData")
load("simpPLPResults.RData")
load("simpWeiResults.RData")
load("WeibullResults.RData")


DICResults <- matrix(0, nrow = 26, ncol = 10)
for(i in c(1, 2, 4, 6, 7, 9, 11:17, 19:21, 25, 26)){
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

DICResults

colSums(DICResults)

expoTot <- 0
for(i in 1:26){
  expoTot <- expoTot + min(DICResults[i,1], DICResults[i,6])
}
expoTot

weiTot <- 0
for(i in 1:26){
  weiTot <- weiTot + min(DICResults[i,2], DICResults[i,7])
}
weiTot

baoTot <- 0
for(i in 1:26){
  baoTot <- baoTot + min(DICResults[i,3], DICResults[i,8])
}
baoTot


#####################################################################

# Posterior Predictive Check

JLTVData <- read.csv("JLTVData.csv")

# Need to load final model draws

JLTVData[, 1] <- as.numeric(JLTVData[, 1])
JLTVData$trun <- JLTVData$MBF == 0

subs <- list()

for(i in 1:20){
  subs[[i]] <- JLTVData[which(JLTVData[,9] == i),]
}

subList <- list()

for(i in 1:20){
  subList[[i]] <- list()
  for(j in 1:8){
    subList[[i]][[j]] <- subs[[i]][which(subs[[i]][,1]==j),c(6,8,3,10)]
  }
}

# getting total miles for each test 
cutoffs <- matrix(0, ncol = 3, nrow = 8)
for(i in 1:8){
  if(i == 5){
    cutoffs[i,1] <- subList[[3]][[i]][1,1]
    cutoffs[i,2] <- subList[[3]][i][[1]][2,1] + subList[[3]][i][[1]][3,1]
    cutoffs[i,3] <- subList[[3]][i][[1]][4,1]
    
  } else {
    cutoffs[i,] <- subList[[3]][i][[1]]$MBF
  }
}

wF.inv <- function(u,lam,beta){
  return((-log(1-u)/lam)^(1/beta))
}

#################################################################################################################
n <- 4000
# choose vehicle type
vNum <- 2

simFailureCounts <- matrix(0, ncol = 3, nrow = n)
sim140Count <- matrix(0, ncol = 3, nrow = n)
trueFailCount <- table(JLTVData[which(JLTVData$Vehicle == vNum),3])
true140Count <- table(JLTVData[which(JLTVData$Vehicle == vNum & JLTVData$MBF < 140),3])

for(i in 1:n){
  # drawNum <- sample(1:1500,1)
  drawNum <- i
  
  for(j in 1:20){
    p1Miles <- 0
    while(p1Miles < cutoffs[vNum,1]){
      newRun <- wF.inv(runif(1), WeibullResults[[j]]$draws[drawNum,vNum + 5], WeibullResults[[j]]$draws[drawNum,3])
      if(newRun < 140){
        sim140Count[i,1] <- sim140Count[i,1] + 1
      }
      p1Miles <- p1Miles + newRun
      simFailureCounts[i,1] <- simFailureCounts[i,1] + 1
    }
    p2Miles <- 0
    while(p2Miles < cutoffs[vNum,2]){
      newRun <- wF.inv(runif(1),
                       WeibullResults[[j]]$draws[drawNum,vNum + 5] * WeibullResults[[j]]$draws[drawNum,4],
                       WeibullResults[[j]]$draws[drawNum,3]) 
      if(newRun < 140){
        sim140Count[i,2] <- sim140Count[i,2] + 1
      }
      p2Miles <- p2Miles + newRun
      simFailureCounts[i,2] <- simFailureCounts[i,2] + 1
    }
    p3Miles <- 0
    while(p3Miles < cutoffs[vNum,3]){
      newRun <- wF.inv(runif(1),
                       WeibullResults[[j]]$draws[drawNum,vNum + 5] * WeibullResults[[j]]$draws[drawNum,4] * WeibullResults[[j]]$draws[drawNum,5],
                       WeibullResults[[j]]$draws[drawNum,3])
      if(newRun < 140){
        sim140Count[i,3] <- sim140Count[i,3] + 1
      }
      p3Miles <- p3Miles + newRun
      simFailureCounts[i,3] <- simFailureCounts[i,3] + 1
    }
  }
}

# total failures 
hist(simFailureCounts[,1], breaks = 50)
abline(v = trueFailCount[1], col = 2)

hist(simFailureCounts[,2], breaks = 50)
abline(v = trueFailCount[2], col = 2)

hist(simFailureCounts[,3], breaks = 50)
abline(v = trueFailCount[3], col = 2)

boxplot(simFailureCounts[,1],simFailureCounts[,2],simFailureCounts[,3],
        main = "Posterior Predictive Simulation Weibull",
        ylab = "Failure Counts",
        xlab = "Phase", 
        names = c("One", "Two", "Three"))
segments(.6, trueFailCount[1], 1.4, trueFailCount[1], col = 2)
segments(1.6, trueFailCount[2], 2.4, trueFailCount[2], col = 2)
segments(2.6, trueFailCount[3], 3.4, trueFailCount[3], col = 2)

# failures in less than 140
hist(sim140Count[,1], breaks = 50)
abline(v = true140Count[1], col = 2)

hist(sim140Count[,2], breaks = 50)
abline(v = true140Count[2], col = 2)

hist(sim140Count[,3], breaks = 50)
abline(v = true140Count[3], col = 2)

boxplot(sim140Count[,1],sim140Count[,2],sim140Count[,3],
        main = "Posterior Predictive Simulation Weibull",
        ylab = "Failure Counts Under 140 Miles",
        xlab = "Phase", 
        names = c("One", "Two", "Three"),
        ylim = c(0, 130))
legend("topright", legend=c("Observed Value"),
       col=c("red"), lty=2, lwd=2)
segments(.6, true140Count[1], 1.4, true140Count[1], col = 2, lwd = 3, lty = 2)
segments(1.6, true140Count[2], 2.4, true140Count[2], col = 2, lwd = 3, lty = 2)
segments(2.6, true140Count[3], 3.4, true140Count[3], col = 2, lwd = 3, lty = 2)


#################################################################################################################

n <- 4000
# choose vehicle type
vNum <- 2

simFailureCounts <- matrix(0, ncol = 3, nrow = n)
sim140Count <- matrix(0, ncol = 3, nrow = n)
trueFailCount <- table(JLTVData[which(JLTVData$Vehicle == vNum),3])
true140Count <- table(JLTVData[which(JLTVData$Vehicle == vNum & JLTVData$MBF < 140),3])

for(i in 1:n){
  # drawNum <- sample(1:1500,1)
  drawNum <- i
  
  for(j in 1:20){
    p1Miles <- 0
    while(p1Miles < cutoffs[vNum,1]){
      newRun <- wF.inv(runif(1), ExpoResults[[j]]$draws[drawNum,vNum + 5], ExpoResults[[j]]$draws[drawNum,3])
      if(newRun < 140){
        sim140Count[i,1] <- sim140Count[i,1] + 1
      }
      p1Miles <- p1Miles + newRun
      simFailureCounts[i,1] <- simFailureCounts[i,1] + 1
    }
    p2Miles <- 0
    while(p2Miles < cutoffs[vNum,2]){
      newRun <- wF.inv(runif(1),
                       ExpoResults[[j]]$draws[drawNum,vNum + 5] * ExpoResults[[j]]$draws[drawNum,4],
                       ExpoResults[[j]]$draws[drawNum,3]) 
      if(newRun < 140){
        sim140Count[i,2] <- sim140Count[i,2] + 1
      }
      p2Miles <- p2Miles + newRun
      simFailureCounts[i,2] <- simFailureCounts[i,2] + 1
    }
    p3Miles <- 0
    while(p3Miles < cutoffs[vNum,3]){
      newRun <- wF.inv(runif(1),
                       ExpoResults[[j]]$draws[drawNum,vNum + 5] * ExpoResults[[j]]$draws[drawNum,4] * ExpoResults[[j]]$draws[drawNum,5],
                       ExpoResults[[j]]$draws[drawNum,3])
      if(newRun < 140){
        sim140Count[i,3] <- sim140Count[i,3] + 1
      }
      p3Miles <- p3Miles + newRun
      simFailureCounts[i,3] <- simFailureCounts[i,3] + 1
    }
  }
}

# total failures 
hist(simFailureCounts[,1], breaks = 50)
abline(v = trueFailCount[1], col = 2)

hist(simFailureCounts[,2], breaks = 50)
abline(v = trueFailCount[2], col = 2)

hist(simFailureCounts[,3], breaks = 50)
abline(v = trueFailCount[3], col = 2)

boxplot(simFailureCounts[,1],simFailureCounts[,2],simFailureCounts[,3],
        main = "Posterior Predictive Simulation Exponential",
        ylab = "Failure Counts",
        xlab = "Phase", 
        names = c("One", "Two", "Three"))
segments(.6, trueFailCount[1], 1.4, trueFailCount[1], col = 2)
segments(1.6, trueFailCount[2], 2.4, trueFailCount[2], col = 2)
segments(2.6, trueFailCount[3], 3.4, trueFailCount[3], col = 2)

# failures in less than 140
hist(sim140Count[,1], breaks = 50)
abline(v = true140Count[1], col = 2)

hist(sim140Count[,2], breaks = 50)
abline(v = true140Count[2], col = 2)

hist(sim140Count[,3], breaks = 50)
abline(v = true140Count[3], col = 2)

boxplot(sim140Count[,1],sim140Count[,2],sim140Count[,3],
        main = "Posterior Predictive Simulation Exponential",
        ylab = "Failure Counts Under 140 Miles",
        xlab = "Phase", 
        names = c("One", "Two", "Three"),
        ylim = c(0, 130))
legend("topright", legend=c("Observed Value"),
       col=c("red"), lty=2, lwd=2)
segments(.6, true140Count[1], 1.4, true140Count[1], col = 2, lwd = 3, lty = 2)
segments(1.6, true140Count[2], 2.4, true140Count[2], col = 2, lwd = 3, lty = 2)
segments(2.6, true140Count[3], 3.4, true140Count[3], col = 2, lwd = 3, lty = 2)



#################################################################################################################


# left out 17 because it made the plot ugly...
plot(density(WeibullResults[[1]]$draws[,3]), ylim = c(0,23), xlim = c(0,0.6), 
     main = "Posterior Plots of Shape Parameter", xlab = "")
for(i in c(2:16, 18:20)){
  lines(density(WeibullResults[[i]]$draws[,3]))
}


