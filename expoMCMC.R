

expoMCMC <- function(data, samples = 5000, shapePriorA = .001,
                             shapePriorB = .001, HyperG1 = .001, HyperG2 = .001,
                             hyperA1 = .001, hyperA2 = .001, hyperT1A = 3,
                             hyperT1B = 3, hyperT2A = 3, hyperT2B = 3,
                             alphaStart = 1, betaStart = 1, theta1Start = 1,
                             theta2Start = 1, tuning = 1,
                             burnin = 1000, thin = 10){
  
  # Save the starting data for imputation
  dataStart <- data
  
  # matrix for keeping MCMC draws for each parameter
  draws <- matrix(0, nrow = samples, ncol = length(data) + 5)
  
  draws[1,1] <- alphaStart
  draws[1,2] <- betaStart
  draws[,3] <- 1
  draws[1,4] <- theta1Start
  draws[1,5] <- theta2Start
  
  # counter for acceptance rate
  acca <- 0
  
  # log posterior function for MCMC draws
  logPost <- function(parm, d) {
    
    lp <- 0
    
    for(k in 6:length(parm)){
      lp <- lp + 
        sum(dexp(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 0 & d[[k-5]]$Phase == 1 & d[[k-5]]$trun == F)], parm[k], log = T)) +
        sum(dexp(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 0 & d[[k-5]]$Phase == 2 & d[[k-5]]$trun == F)], parm[k]*parm[4], log = T)) +
        sum(dexp(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 0 & d[[k-5]]$Phase == 3 & d[[k-5]]$trun == F)], parm[k]*parm[4]*parm[5], log = T)) -
        sum(parm[k]*(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 1 & d[[k-5]]$Phase == 1 & d[[k-5]]$trun == F)])) -
        sum(parm[k]*parm[4]*(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 1 & d[[k-5]]$Phase == 2 & d[[k-5]]$trun == F)])) -
        sum(parm[k]*parm[4]*parm[5]*(d[[k-5]]$MBF[which(d[[k-5]]$Censor == 1 & d[[k-5]]$Phase == 3 & d[[k-5]]$trun == F)])) +
        sum(parm[k]*(d[[k-5]]$Lower[which(d[[k-5]]$Phase == 1 & d[[k-5]]$trun == T)]) - 
              parm[k]*(d[[k-5]]$Upper[which(d[[k-5]]$Phase == 1 & d[[k-5]]$trun == T)])) +
        sum(parm[k]*parm[4]*(d[[k-5]]$Lower[which(d[[k-5]]$Phase == 2 & d[[k-5]]$trun == T)]) - 
              parm[k]*parm[4]*(d[[k-5]]$Upper[which(d[[k-5]]$Phase == 2 & d[[k-5]]$trun == T)])) + 
        sum(parm[k]*parm[4]*parm[5]*(d[[k-5]]$Lower[which(d[[k-5]]$Phase == 3 & d[[k-5]]$trun == T)]) - 
              parm[k]*parm[4]*parm[5]*(d[[k-5]]$Upper[which(d[[k-5]]$Phase == 3 & d[[k-5]]$trun == T)])) +
        dgamma(parm[k], parm[1], parm[2], log=T)
    }
    
    lp <- lp + 
      dgamma(parm[1], hyperA1, hyperA2, log=T) + 
      dgamma(parm[2],HyperG1, HyperG2, log=T) +
      dgamma(parm[3], shapePriorA, shapePriorB, log=T) +
      dgamma(parm[4],hyperT1A, hyperT1B, log=T) +
      dgamma(parm[5],hyperT2A, hyperT2B, log=T)
    
    return(lp)
  }
  
  # phase obs counts 
  phase2Count <- 0
  phase3Count <- 0
  for(i in 1:length(data)){
    phase2Count <- phase2Count + length(which(data[[i]]$Phase == 2 & data[[i]]$Censor == 0))
    phase3Count <- phase3Count + length(which(data[[i]]$Phase == 3 & data[[i]]$Censor == 0))
  }
  
  # MCMC draws
  for (i in 2:samples) {
    
    for(j in 6:ncol(draws)){
      draws[i,j] <- rgamma(1,sum(data[[j-5]]$Censor == 0)+draws[i-1,1],
                           sum(data[[j-5]]$MBF[which(data[[j-5]]$Phase == 1)],
                               (draws[i-1,4] * data[[j-5]]$MBF[which(data[[j-5]]$Phase == 2)]),
                               (draws[i-1,4] * draws[i-1,5] * data[[j-5]]$MBF[which(data[[j-5]]$Phase == 3)]))
                           + draws[i-1,2])
    }
    
    draws[i,2] <- rgamma(1, length(data)*draws[i-1,1] + HyperG1,
                         sum(draws[i,6:ncol(draws)]) + HyperG2)
    
    # Phase sums
    phase2Sum <- 0
    phase3Sum4 <- 0
    phase3Sum5 <- 0
    for(k in 1:length(data)){
      phase2Sum <- phase2Sum + draws[i,k+5] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 2)])
      
      # for the first theta
      phase3Sum4 <- phase3Sum4 + draws[i,k+5] * draws[i-1,5] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 3)])
      
      # for second theta
      phase3Sum5 <- phase3Sum5 + draws[i,k+5] * draws[i-1,4] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 3)])
    }
    
    draws[i,4] <- rgamma(1, phase2Count + phase3Count + hyperT1A, phase2Sum + phase3Sum4 + hyperT1B)
    draws[i,5] <- rgamma(1, phase3Count + hyperT2A, phase3Sum5 + hyperT2B)
    
    
    ################################
    # metropolis hastings
    ################################
    
    # Auto-adjusting Tuning Params
    if((i %% 500) == 0){
      if(acca/i > .55){
        if(acca/i > .7){
          tuning <- tuning*2
        } else{
          tuning <- tuning*1.5
        }
      } else if (acca/i < .3) {
        if(acca/i < .2){
          tuning <- tuning/2
        } else {
          tuning <- tuning/1.5
        }
      }
    }
    
    # Sample from alpha 
    draws[i,1] <- draws[i-1,1]
    draws[i,3] <- draws[i-1,3]
    astar <- rnorm(1, draws[i-1,1], tuning)
    if (astar > 0) {
      lnew <- logPost(c(astar, draws[i,2:ncol(draws)]), data)
      lold <- logPost(c(draws[i-1,1], draws[i,2:ncol(draws)]), data)
      if(is.finite(lnew - lold)){
        if (lnew - lold > log(runif(1))) {
          draws[i,1] <- astar
          acca <- acca + 1
        }
      }
    }
  }
  
  finalDraws <- draws[seq(from = burnin + 1, to = samples, by = thin),]
  
  # DIC
  d <- rep(0, nrow(finalDraws))
  for(i in 1:nrow(finalDraws)){
    for(j in 6:ncol(finalDraws)){
      d[i] <- d[i] - 
        2*(sum(dexp(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 1 & data[[j-5]]$trun == F)], finalDraws[i,j], log = T)) +
             sum(dexp(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 2 & data[[j-5]]$trun == F)], finalDraws[i,j]*finalDraws[i,4], log = T)) +
             sum(dexp(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 3 & data[[j-5]]$trun == F)], finalDraws[i,j]*finalDraws[i,4]*finalDraws[i,5], log = T)) -
             sum(finalDraws[i,j]*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 1 & data[[j-5]]$trun == F)])) -
             sum(finalDraws[i,j]*finalDraws[i,4]*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 2 & data[[j-5]]$trun == F)])) -
             sum(finalDraws[i,j]*finalDraws[i,4]*finalDraws[i,5]*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 3 & data[[j-5]]$trun == F)])) +
             sum(finalDraws[i,j]*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 1 & data[[j-5]]$trun == T)]) - 
                  finalDraws[i,j]*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 1 & data[[j-5]]$trun == T)])) +
             sum(finalDraws[i,j]*finalDraws[i,4]*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 2 & data[[j-5]]$trun == T)]) - 
                  finalDraws[i,j]*finalDraws[i,4]*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 2 & data[[j-5]]$trun == T)])) + 
             sum(finalDraws[i,j]*finalDraws[i,4]*finalDraws[i,5]*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 3 & data[[j-5]]$trun == T)]) - 
                  finalDraws[i,j]*finalDraws[i,4]*finalDraws[i,5]*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 3 & data[[j-5]]$trun == T)])))
    }
  }
    
  davg <- mean(d)
  dthetahat <- 0
  for(j in 6:ncol(finalDraws)){
    dthetahat <- dthetahat - 2*(sum(dexp(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 1 & data[[j-5]]$trun == F)], mean(finalDraws[,j]), log = T)) +
                                  sum(dexp(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 2 & data[[j-5]]$trun == F)], mean(finalDraws[,j])*mean(finalDraws[,4]), log = T)) +
                                  sum(dexp(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 0 & data[[j-5]]$Phase == 3 & data[[j-5]]$trun == F)], mean(finalDraws[,j])*mean(finalDraws[,4])*mean(finalDraws[,5]), log = T)) -
                                  sum(mean(finalDraws[,j])*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 1 & data[[j-5]]$trun == F)])) -
                                  sum(mean(finalDraws[,j])*mean(finalDraws[,4])*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 2 & data[[j-5]]$trun == F)])) -
                                  sum(mean(finalDraws[,j])*mean(finalDraws[,4])*mean(finalDraws[,5])*(data[[j-5]]$MBF[which(data[[j-5]]$Censor == 1 & data[[j-5]]$Phase == 3 & data[[j-5]]$trun == F)])) +
                                  sum(mean(finalDraws[,j])*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 1 & data[[j-5]]$trun == T)]) - 
                                      mean(finalDraws[,j])*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 1 & data[[j-5]]$trun == T)])) +
                                  sum(mean(finalDraws[,j])*mean(finalDraws[,4])*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 2 & data[[j-5]]$trun == T)]) - 
                                      mean(finalDraws[,j])*mean(finalDraws[,4])*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 2 & data[[j-5]]$trun == T)])) + 
                                  sum(mean(finalDraws[,j])*mean(finalDraws[,4])*mean(finalDraws[,5])*(data[[j-5]]$Lower[which(data[[j-5]]$Phase == 3 & data[[j-5]]$trun == T)]) - 
                                      mean(finalDraws[,j])*mean(finalDraws[,4])*mean(finalDraws[,5])*(data[[j-5]]$Upper[which(data[[j-5]]$Phase == 3 & data[[j-5]]$trun == T)])))
  }
    
  pd <- davg - dthetahat
  dic <- davg + pd
  
  return(list(draws = finalDraws,
              acceptance = acca/samples, 
              DIC = dic,
              PD = pd,
              Deviance = d))
}