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
    if(all(subList[[i]][[j]]$trun == F)){
    } else {
      for(k in 1:nrow(subList[[i]][[j]]))
        if(subList[[i]][[j]][k,4] == T){
          if(k == 1){
          } else {
            subList[[i]][[j]][k,5] <- subList[[i]][[j]][k,7]
            if(subList[[i]][[j]][k-1,2] == 1 & subList[[i]][[j]][k-1,3] != 3){
              subList[[i]][[j]][k,3] <- subList[[i]][[j]][k-1,3] + 1
            }
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
        }
    }
}


