
simpExpoMCMC <- function(data, samples = 50000, shapePriorA = .001,
                     shapePriorB = .001, lamPriorA = .001, lamPriorB = .001,
                     theta1PriorA = .001, theta1PriorB = .001, theta2PriorA = .001,
                     theta2PriorB = .001, theta1Start = 1, theta2Start = 1, tuning = 1,
                     burnin = 35000, thin = 10){
  
  # matrix for keeping MCMC draws for each parameter
  lam_draws <- rep(0, samples)
  lam_draws[1] <- 1
  
  theta1_draws <- rep(0, samples)
  theta1_draws[1] <- theta1Start
  
  theta2_draws <- rep(0, samples)
  theta2_draws[1] <- theta2Start
  
  # phase obs counts 
  phase2Count <- 0
  phase3Count <- 0
  for(i in 1:length(data)){
    phase2Count <- phase2Count + length(which(data[[i]]$Phase == 2 & data[[i]]$Censor == 0))
    phase3Count <- phase3Count + length(which(data[[i]]$Phase == 3 & data[[i]]$Censor == 0))
  }
  
  noCenCount <- 0
  for(i in 1:length(data)){
    noCenCount <- noCenCount + sum(data[[i]]$Censor == 0)
  }
  
  # MCMC draws
  for (i in 2:samples) {
    
    datSum <- 0
    for(j in 1:length(data)){
      datSum <- datSum + sum(data[[j]]$MBF[which(data[[j]]$Phase == 1)],
                             (theta1_draws[i-1] * data[[j]]$MBF[which(data[[j]]$Phase == 2)]),
                             (theta1_draws[i-1] * theta2_draws[i-1] * data[[j]]$MBF[which(data[[j]]$Phase == 3)]))
    }
    
    lam_draws[i] <- rgamma(1, noCenCount + lamPriorA, datSum + lamPriorB)

    # Phase sums
    phase2Sum <- 0
    phase3Sum4 <- 0
    phase3Sum5 <- 0
    for(k in 1:length(data)){
      phase2Sum <- phase2Sum + lam_draws[i] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 2)])
      
      # for the first theta
      phase3Sum4 <- phase3Sum4 + theta2_draws[i-1] * lam_draws[i] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 3)])
      
      # for second theta
      phase3Sum5 <- phase3Sum5 + theta1_draws[i-1] * lam_draws[i] * sum(data[[k]]$MBF[which(data[[k]]$Phase == 3)])
    }
    
    theta1_draws[i] <- rgamma(1, phase2Count + phase3Count + theta1PriorA, phase2Sum + phase3Sum4 + theta1PriorB)
    theta2_draws[i] <- rgamma(1, phase3Count + theta2PriorA, phase3Sum5 + theta2PriorB)
  }
  
  lam_finalDraws <- lam_draws[seq(from = burnin + 1, to = samples, by = thin)]
  theta1_finalDraws <- theta1_draws[seq(from = burnin + 1, to = samples, by = thin)]
  theta2_finalDraws <- theta2_draws[seq(from = burnin + 1, to = samples, by = thin)]
  
  # DIC
  d <- rep(0, length(lam_finalDraws))
  for(i in 1:length(lam_finalDraws)){
    for(j in 1:length(data)){
      d[i] <- d[i] - 
        2*(sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1 & data[[j]]$trun == F)], lam_finalDraws[i], log = T)) +
             sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2 & data[[j]]$trun == F)], lam_finalDraws[i]*theta1_finalDraws[i], log = T)) +
             sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3 & data[[j]]$trun == F)], lam_finalDraws[i]*theta1_finalDraws[i]*theta2_finalDraws[i], log = T)) -
             sum(lam_finalDraws[i]*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 1 & data[[j]]$trun == F)])) -
             sum(lam_finalDraws[i]*theta1_finalDraws[i]*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 2 & data[[j]]$trun == F)])) -
             sum(lam_finalDraws[i]*theta1_finalDraws[i]*theta2_finalDraws[i]*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 3 & data[[j]]$trun == F)])) +
             sum(lam_finalDraws[i]*(data[[j]]$Lower[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)]) - 
                   lam_finalDraws[i]*(data[[j]]$Upper[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)])) +
             sum(lam_finalDraws[i]*theta1_finalDraws[i]*(data[[j]]$Lower[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)]) - 
                   lam_finalDraws[i]*theta1_finalDraws[i]*(data[[j]]$Upper[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)])) + 
             sum(lam_finalDraws[i]*theta1_finalDraws[i]*theta2_finalDraws[i]*(data[[j]]$Lower[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)]) - 
                   lam_finalDraws[i]*theta1_finalDraws[i]*theta2_finalDraws[i]*(data[[j]]$Upper[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)])))
    }
  }
  
  davg <- mean(d)
  dthetahat <- 0
  for(j in 1:length(data)){
    dthetahat <- dthetahat - 2*(sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1 & data[[j]]$trun == F)], mean(lam_finalDraws), log = T)) +
                                  sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2 & data[[j]]$trun == F)], mean(lam_finalDraws)*mean(theta1_finalDraws), log = T)) +
                                  sum(dexp(data[[j]]$MBF[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3 & data[[j]]$trun == F)], mean(lam_finalDraws)*mean(theta1_finalDraws)*mean(theta2_finalDraws), log = T)) -
                                  sum(mean(lam_finalDraws)*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 1 & data[[j]]$trun == F)])) -
                                  sum(mean(lam_finalDraws)*mean(theta1_finalDraws)*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 2 & data[[j]]$trun == F)])) -
                                  sum(mean(lam_finalDraws)*mean(theta1_finalDraws)*mean(theta2_finalDraws)*(data[[j]]$MBF[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 3 & data[[j]]$trun == F)])) +
                                  sum(mean(lam_finalDraws)*(data[[j]]$Lower[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)]) - 
                                        mean(lam_finalDraws)*(data[[j]]$Upper[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)])) +
                                  sum(mean(lam_finalDraws)*mean(theta1_finalDraws)*(data[[j]]$Lower[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)]) - 
                                        mean(lam_finalDraws)*mean(theta1_finalDraws)*(data[[j]]$Upper[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)])) + 
                                  sum(mean(lam_finalDraws)*mean(theta1_finalDraws)*mean(theta2_finalDraws)*(data[[j]]$Lower[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)]) - 
                                        mean(lam_finalDraws)*mean(theta1_finalDraws)*mean(theta2_finalDraws)*(data[[j]]$Upper[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)])))
  }
  
  pd <- davg - dthetahat
  dic <- davg + pd
  
  return(list(lam_draws = lam_finalDraws,
              theta1_draws = theta1_finalDraws,
              theta2_draws = theta2_finalDraws,
              DIC = dic,
              PD = pd,
              Deviance = d))
}