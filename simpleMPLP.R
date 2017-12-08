logMPLPsum <- function(U, L, kappa, shape, scale){
  d <- 0
  if(length(U) > 0){
    d <- sum(-log(gamma(kappa)) + log(shape) + kappa*log(scale) + (shape - 1)*log(U) + (kappa - 1)*(log((U^shape - L^shape))))
  }
  return(d)
}

logMPLPInterval <- function(L, U, kappa, shape, scale){
  d <- 0
  if(length(U) > 0){
    d <- sum(-log(gamma(kappa)) + (kappa)*(log(scale*(U^shape - L^shape))) - scale*(U^shape - L^shape))
  } 
  return(d)
}


simpMPLPMCMC <- function(data, samples = 100000, shapePriorA = .001,
                     shapePriorB = .001, priorT1A = 3, priorT1B = 3, priorT2A = 3, priorT2B = 3,
                     theta1Start = 1, theta2Start = 1, shapeStart = 1, kappaStart = 1, tuningK = 1,
                     tuningS = 1, burnin = 75000, thin = 10, lamPriorA = .001, lamPriorB = .001, 
                     priorKA = .001, priorKB = .001){
  
  # matrix for keeping MCMC draws for each parameter
  draws <- matrix(0, nrow = samples, ncol = 7)
  
  draws[1,3] <- shapeStart
  draws[1,4] <- theta1Start
  draws[1,5] <- theta2Start
  draws[1,6] <- kappaStart
  
  # counter for acceptance rate
  accs <- 0
  acck <- 0
  
  # log posterior function for MCMC draws
  logPost <- function(parm, d) {
    
    lp <- 0
    
    for(k in 1:length(d)){
      
      lp <- lp + 
        logMPLPsum(d[[k]]$Upper[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 1 & d[[k]]$trun == F)], d[[k]]$Lower[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 1 & d[[k]]$trun == F)], parm[6], parm[3], parm[7]) +
        logMPLPsum(d[[k]]$Upper[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 2 & d[[k]]$trun == F)], d[[k]]$Lower[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 2 & d[[k]]$trun == F)], parm[6], parm[3], parm[7]*parm[4]) +
        logMPLPsum(d[[k]]$Upper[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 3 & d[[k]]$trun == F)], d[[k]]$Lower[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 3 & d[[k]]$trun == F)], parm[6], parm[3], parm[7]*parm[4]*parm[5]) -
        sum(parm[7]*(d[[k]]$totalMiles[which(d[[k]]$Censor == 1 & d[[k]]$Phase == 1 & d[[k]]$trun == F)]^parm[3])) -
        sum(parm[7]*parm[4]*(d[[k]]$totalMiles[which(d[[k]]$Censor == 1 & d[[k]]$Phase == 2 & d[[k]]$trun == F)]^parm[3])) -
        sum(parm[7]*parm[4]*parm[5]*(d[[k]]$totalMiles[which(d[[k]]$Censor == 1 & d[[k]]$Phase == 3 & d[[k]]$trun == F)]^parm[3])) +
        logMPLPInterval(d[[k]]$Lower[which(d[[k]]$Phase == 1 & d[[k]]$trun == T)], d[[k]]$Upper[which(d[[k]]$Phase == 1 & d[[k]]$trun == T)], parm[6], parm[3], parm[7]) +
        logMPLPInterval(d[[k]]$Lower[which(d[[k]]$Phase == 2 & d[[k]]$trun == T)], d[[k]]$Upper[which(d[[k]]$Phase == 2 & d[[k]]$trun == T)], parm[6], parm[3], parm[7]*parm[4]) +
        logMPLPInterval(d[[k]]$Lower[which(d[[k]]$Phase == 3 & d[[k]]$trun == T)], d[[k]]$Upper[which(d[[k]]$Phase == 3 & d[[k]]$trun == T)], parm[6], parm[3], parm[7]*parm[4]*parm[5]) 
    }
    lp <- lp + 
      dgamma(parm[7],lamPriorA, lamPriorB, log=T) +
      dgamma(parm[3], shapePriorA, shapePriorB, log=T) +
      dgamma(parm[4],priorT1A, priorT1B, log=T) +
      dgamma(parm[5],priorT2A, priorT2B, log=T) + 
      dgamma(parm[6],priorKA, priorKB, log=T)
    
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
    
    noCenCount <- 0
    for(j in 1:length(data)){
      noCenCount <- noCenCount + sum(data[[j]]$Censor == 0)*draws[i-1,6]
    }

    datSum <- 0
    for(j in 1:length(data)){
      datSum <- datSum + sum(data[[j]]$totalMiles[which(data[[j]]$Phase == 1 & data[[j]]$Censor == 1 & data[[j]]$trun == F)]^draws[i-1,3],
                             (draws[i-1,4] * data[[j]]$totalMiles[which(data[[j]]$Phase == 2 & data[[j]]$Censor == 1 & data[[j]]$trun == F)]^draws[i-1,3]),
                             (draws[i-1,4] * draws[i-1,5] * data[[j]]$totalMiles[which(data[[j]]$Phase == 3 & data[[j]]$Censor == 1 & data[[j]]$trun == F)]^draws[i-1,3])) +
        sum(data[[j]]$Upper[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)]^draws[i-1,3] - 
              data[[j]]$Lower[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)]^draws[i-1,3],
            draws[i-1,4] * (data[[j]]$Upper[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)]^draws[i-1,3] - 
                              data[[j]]$Lower[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)]^draws[i-1,3]),
            draws[i-1,4] * draws[i-1,5] *(data[[j]]$Upper[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)]^draws[i-1,3] - 
                                            data[[j]]$Lower[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)]^draws[i-1,3]))
    }
    
    draws[i,7] <- rgamma(1, noCenCount + lamPriorA, datSum + lamPriorB)
    
    # Phase sums
    phase2Sum <- 0
    phase3Sum4 <- 0
    phase3Sum5 <- 0
    for(k in 1:length(data)){
      phase2Sum <- phase2Sum + draws[i,7] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 2 & data[[k]]$Censor == 1 & data[[k]]$trun == F)]^draws[i-1,3],
                                                  data[[k]]$Upper[which(data[[k]]$Phase == 2 & data[[k]]$trun == T)]^draws[i-1,3] -
                                                    data[[k]]$Lower[which(data[[k]]$Phase == 2 & data[[k]]$trun == T)]^draws[i-1,3])
      
      # for the first theta
      phase3Sum4 <- phase3Sum4 + draws[i,7] * draws[i-1,5] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 3 & data[[k]]$Censor == 1 & data[[k]]$trun == F)]^draws[i-1,3],
                                                                   data[[k]]$Upper[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3] -
                                                                     data[[k]]$Lower[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3])
      
      # for second theta
      phase3Sum5 <- phase3Sum5 + draws[i,7] * draws[i-1,4] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 3 & data[[k]]$Censor == 1 & data[[k]]$trun == F)]^draws[i-1,3],
                                                                   data[[k]]$Upper[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3] -
                                                                     data[[k]]$Lower[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3])
    }
    
    draws[i,4] <- rgamma(1, (phase2Count + phase3Count)*draws[i-1,6] + priorT1A, phase2Sum + phase3Sum4 + priorT1B)
    draws[i,5] <- rgamma(1, phase3Count*draws[i-1,6] + priorT2A, phase3Sum5 + priorT2B)
    
    
    ################################
    # metropolis hastings
    ################################
    
    # Auto-adjusting Tuning Params
    if((i %% 500) == 0){
      if(accs/i > .55){
        if(accs/i > .7){
          tuningS <- tuningS*2
        } else {
          tuningS <- tuningS*1.5
        }
      } else if (accs/i < .3) { 
        if (accs/i < .2){
          tuningS <- tuningS/2
        } else {
          tuningS <- tuningS/1.5
        }
      }
      if(acck/i > .55){
        if(acck/i > .7){
          tuningK <- tuningK*2
        } else {
          tuningK <- tuningK*1.5
        }
      } else if (acck/i < .3) { 
        if (acck/i < .2){
          tuningK <- tuningK/2
        } else {
          tuningK <- tuningK/1.5
        }
      }
    }
    
    # Sample from alpha 
    draws[i,3] <- draws[i-1,3]
    draws[i,6] <- draws[i-1,6]
    # Sample from shape
    sstar <- rnorm(1, draws[i-1,3], tuningS)
    if (sstar > 0) {
      lnew <- logPost(c( draws[i,1:2], sstar, draws[i,4:ncol(draws)]), data)
      lold <- logPost(c( draws[i,1:2], draws[i-1,3], draws[i,4:ncol(draws)]), data)
      if(is.finite(lnew - lold)){
        if (lnew - lold > log(runif(1))) {
          draws[i,3] <- sstar
          accs <- accs + 1
        }
      }
    }
    
    # Sample from kappa
    kstar <- rnorm(1, draws[i-1,6], tuningK)
    if (kstar > 0) {
      lnew <- logPost(c( draws[i,1:5], kstar, draws[i,7]), data)
      lold <- logPost(c( draws[i,1:5], draws[i-1,6], draws[i,7]), data)
      if(is.finite(lnew - lold)){
        if (lnew - lold > log(runif(1))) {
          draws[i,6] <- kstar
          acck <- acck + 1
        }
      }
    }
  }
  
  finalDraws <- draws[seq(from = burnin + 1, to = samples, by = thin),]
  
  # DIC
  d <- rep(0, nrow(finalDraws))
  for(i in 1:nrow(finalDraws)){
    for(j in 1:length(data)){
      d[i] <- d[i] - 
        2*(logMPLPsum(data[[j]]$Upper[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1 & data[[j]]$trun == F)], data[[j]]$Lower[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1 & data[[j]]$trun == F)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,7]) +
             logMPLPsum(data[[j]]$Upper[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2 & data[[j]]$trun == F)], data[[j]]$Lower[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2 & data[[j]]$trun == F)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,7]*finalDraws[i,4]) +
             logMPLPsum(data[[j]]$Upper[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3 & data[[j]]$trun == F)], data[[j]]$Lower[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3 & data[[j]]$trun == F)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,7]*finalDraws[i,4]*finalDraws[i,5]) -
             sum(finalDraws[i,7]*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 1 & data[[j]]$trun == F)]^finalDraws[i,3])) -
             sum(finalDraws[i,7]*finalDraws[i,4]*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 2 & data[[j]]$trun == F)]^finalDraws[i,3])) -
             sum(finalDraws[i,7]*finalDraws[i,4]*finalDraws[i,5]*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 3 & data[[j]]$trun == F)]^finalDraws[i,3])) +
             logMPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,7]) +
             logMPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,7]*finalDraws[i,4]) +
             logMPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,7]*finalDraws[i,4]*finalDraws[i,5]))
    }
  }
  
  davg <- mean(d)
  dthetahat <- 0
  for(j in 1:length(data)){
    dthetahat <- dthetahat - 2*(logMPLPsum(data[[j]]$Upper[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1 & data[[j]]$trun == F)], data[[j]]$Lower[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1 & data[[j]]$trun == F)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,7])) +
                                  logMPLPsum(data[[j]]$Upper[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2 & data[[j]]$trun == F)], data[[j]]$Lower[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2 & data[[j]]$trun == F)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,7])*mean(finalDraws[,4])) +
                                  logMPLPsum(data[[j]]$Upper[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3 & data[[j]]$trun == F)], data[[j]]$Lower[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3 & data[[j]]$trun == F)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,7])*mean(finalDraws[,4])*mean(finalDraws[,5])) -
                                  sum(mean(finalDraws[,7])*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 1 & data[[j]]$trun == F)]^mean(finalDraws[,3]))) -
                                  sum(mean(finalDraws[,7])*mean(finalDraws[,4])*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 2 & data[[j]]$trun == F)]^mean(finalDraws[,3]))) -
                                  sum(mean(finalDraws[,7])*mean(finalDraws[,4])*mean(finalDraws[,5])*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 3 & data[[j]]$trun == F)]^mean(finalDraws[,3]))) +
                                  logMPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,7])) +
                                  logMPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,7])*mean(finalDraws[,4])) +
                                  logMPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,7])*mean(finalDraws[,4])*mean(finalDraws[,5])))
  }
  
  pd <- davg - dthetahat
  dic <- davg + pd
  
  return(list(draws = finalDraws,
              acceptanceShape = accs/samples,
              acceptanceKappa = acck/samples,
              DIC = dic,
              PD = pd,
              Deviance = d))
}