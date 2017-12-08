
logMPLPsum <- function(U, L, kappa, shape, scale){
  d <- 0
  if(length(U) > 0){
    d <- sum(-log(gamma(kappa)) + log(shape) + kappa*log(scale) + (shape - 1)*log(U) + (kappa - 1)*(log((U^shape - L^shape))))
    return(d)
  } else {
    return(0)
  }

}

logMPLPInterval <- function(L, U, kappa, shape, scale){
  d <- 0
  if(length(U) > 0){
    d <- sum(-log(gamma(kappa)) + (kappa)*(log(scale*(U^shape - L^shape))) - scale*(U^shape - L^shape))
    return(d)
  } else {
    return(0)
  }

}


MPLPMCMC <- function(data, samples = 100000, shapePriorA = .001,
                    shapePriorB = .001, HyperG1 = .001, HyperG2 = .001,
                    hyperA1 = .001, hyperA2 = .001, hyperT1A = 3,
                    hyperT1B = 3, hyperT2A = 3, hyperT2B = 3, hyperKA = .001, hyperKB = .001,
                    alphaStart = 1, betaStart = 1, theta1Start = 1,
                    theta2Start = 1, shapeStart = 1, kappaStart = 1, tuningA = 1,
                    tuningS = 1, tuningK = 1, burnin = 75000, thin = 10){
  
  # Save the starting data for imputation
  dataStart <- data
  
  # matrix for keeping MCMC draws for each parameter
  draws <- matrix(0, nrow = samples, ncol = length(data) + 6)
  
  draws[1,1] <- alphaStart
  draws[1,2] <- betaStart
  draws[1,3] <- shapeStart
  draws[1,4] <- theta1Start
  draws[1,5] <- theta2Start
  draws[1,6] <- kappaStart
  
  # counter for acceptance rate
  acca <- 0
  accs <- 0
  acck <- 0
  
  # log posterior function for MCMC draws
  logPost <- function(parm, d) {

    lp <- 0
    
    for(k in 7:length(parm)){
      
      lp <- lp + 
        logMPLPsum(d[[k-6]]$Upper[which(d[[k-6]]$Censor == 0 & d[[k-6]]$Phase == 1 & d[[k-6]]$trun == F)], d[[k-6]]$Lower[which(d[[k-6]]$Censor == 0 & d[[k-6]]$Phase == 1 & d[[k-6]]$trun == F)], parm[6], parm[3], parm[k]) +
        logMPLPsum(d[[k-6]]$Upper[which(d[[k-6]]$Censor == 0 & d[[k-6]]$Phase == 2 & d[[k-6]]$trun == F)], d[[k-6]]$Lower[which(d[[k-6]]$Censor == 0 & d[[k-6]]$Phase == 2 & d[[k-6]]$trun == F)], parm[6], parm[3], parm[k]*parm[4]) +
        logMPLPsum(d[[k-6]]$Upper[which(d[[k-6]]$Censor == 0 & d[[k-6]]$Phase == 3 & d[[k-6]]$trun == F)], d[[k-6]]$Lower[which(d[[k-6]]$Censor == 0 & d[[k-6]]$Phase == 3 & d[[k-6]]$trun == F)], parm[6], parm[3], parm[k]*parm[4]*parm[5]) -
        sum(parm[k]*(d[[k-6]]$totalMiles[which(d[[k-6]]$Censor == 1 & d[[k-6]]$Phase == 1 & d[[k-6]]$trun == F)]^parm[3])) -
        sum(parm[k]*parm[4]*(d[[k-6]]$totalMiles[which(d[[k-6]]$Censor == 1 & d[[k-6]]$Phase == 2 & d[[k-6]]$trun == F)]^parm[3])) -
        sum(parm[k]*parm[4]*parm[5]*(d[[k-6]]$totalMiles[which(d[[k-6]]$Censor == 1 & d[[k-6]]$Phase == 3 & d[[k-6]]$trun == F)]^parm[3])) +
        logMPLPInterval(d[[k-6]]$Lower[which(d[[k-6]]$Phase == 1 & d[[k-6]]$trun == T)], d[[k-6]]$Upper[which(d[[k-6]]$Phase == 1 & d[[k-6]]$trun == T)], parm[6], parm[3], parm[k]) +
        logMPLPInterval(d[[k-6]]$Lower[which(d[[k-6]]$Phase == 2 & d[[k-6]]$trun == T)], d[[k-6]]$Upper[which(d[[k-6]]$Phase == 2 & d[[k-6]]$trun == T)], parm[6], parm[3], parm[k]*parm[4]) +
        logMPLPInterval(d[[k-6]]$Lower[which(d[[k-6]]$Phase == 3 & d[[k-6]]$trun == T)], d[[k-6]]$Upper[which(d[[k-6]]$Phase == 3 & d[[k-6]]$trun == T)], parm[6], parm[3], parm[k]*parm[4]*parm[5]) +
        dgamma(parm[k], parm[1], parm[2], log=T)
    }
    lp <- lp + 
      dgamma(parm[1], hyperA1, hyperA2, log=T) + 
      dgamma(parm[2],HyperG1, HyperG2, log=T) +
      dgamma(parm[3], shapePriorA, shapePriorB, log=T) +
      dgamma(parm[4],hyperT1A, hyperT1B, log=T) +
      dgamma(parm[5],hyperT2A, hyperT2B, log=T) + 
      dgamma(parm[6],hyperKA, hyperKB, log=T)

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
    
    for(j in 7:ncol(draws)){
      draws[i,j] <- rgamma(1,sum(data[[j-6]]$Censor == 0)*draws[i-1,6] + draws[i-1,1],
                           sum(data[[j-6]]$totalMiles[which(data[[j-6]]$Phase == 1 & data[[j-6]]$Censor == 1 & data[[j-6]]$trun == F)]^draws[i-1,3],
                               (draws[i-1,4] * data[[j-6]]$totalMiles[which(data[[j-6]]$Phase == 2 & data[[j-6]]$Censor == 1 & data[[j-6]]$trun == F)]^draws[i-1,3]),
                               (draws[i-1,4] * draws[i-1,5] * data[[j-6]]$totalMiles[which(data[[j-6]]$Phase == 3 & data[[j-6]]$Censor == 1 & data[[j-6]]$trun == F)]^draws[i-1,3])) +
                             sum(data[[j-6]]$Upper[which(data[[j-6]]$Phase == 1 & data[[j-6]]$trun == T)]^draws[i-1,3] - 
                                   data[[j-6]]$Lower[which(data[[j-6]]$Phase == 1 & data[[j-6]]$trun == T)]^draws[i-1,3],
                                 draws[i-1,4] * (data[[j-6]]$Upper[which(data[[j-6]]$Phase == 2 & data[[j-6]]$trun == T)]^draws[i-1,3] - 
                                                   data[[j-6]]$Lower[which(data[[j-6]]$Phase == 2 & data[[j-6]]$trun == T)]^draws[i-1,3]),
                                 draws[i-1,4] * draws[i-1,5] *(data[[j-6]]$Upper[which(data[[j-6]]$Phase == 3 & data[[j-6]]$trun == T)]^draws[i-1,3] - 
                                                                 data[[j-6]]$Lower[which(data[[j-6]]$Phase == 3 & data[[j-6]]$trun == T)]^draws[i-1,3]))  
                           + draws[i-1,2])
    }
    
    draws[i,2] <- rgamma(1, length(data)*draws[i-1,1] + HyperG1,
                         sum(draws[i,6:ncol(draws)]) + HyperG2)
    
    
    
    # Phase sums
    phase2Sum <- 0
    phase3Sum4 <- 0
    phase3Sum5 <- 0
    for(k in 1:length(data)){
      phase2Sum <- phase2Sum + draws[i,k+6] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 2 & data[[k]]$Censor == 1 & data[[k]]$trun == F)]^draws[i-1,3],
                                                  data[[k]]$Upper[which(data[[k]]$Phase == 2 & data[[k]]$trun == T)]^draws[i-1,3] -
                                                    data[[k]]$Lower[which(data[[k]]$Phase == 2 & data[[k]]$trun == T)]^draws[i-1,3])
      
      # for the first theta
      phase3Sum4 <- phase3Sum4 + draws[i,k+6] * draws[i-1,5] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 3 & data[[k]]$Censor == 1 & data[[k]]$trun == F)]^draws[i-1,3],
                                                                   data[[k]]$Upper[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3] -
                                                                     data[[k]]$Lower[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3])
      
      # for second theta
      phase3Sum5 <- phase3Sum5 + draws[i,k+6] * draws[i-1,4] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 3 & data[[k]]$Censor == 1 & data[[k]]$trun == F)]^draws[i-1,3],
                                                                   data[[k]]$Upper[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3] -
                                                                     data[[k]]$Lower[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3])
    }
    
    draws[i,4] <- rgamma(1, (phase2Count + phase3Count)*draws[i-1,6] + hyperT1A, phase2Sum + phase3Sum4 + hyperT1B)
    draws[i,5] <- rgamma(1, phase3Count*draws[i-1,6] + hyperT2A, phase3Sum5 + hyperT2B)
    
    
    ################################
    # metropolis hastings
    ################################
    
    # Auto-adjusting Tuning Params
    if((i %% 500) == 0){
      if(acca/i > .55){
        if(acca/i > .7){
          tuningA <- tuningA*2
        } else{
          tuningA <- tuningA*1.5
        }
      } else if (acca/i < .3) {
        if(acca/i < .2){
          tuningA <- tuningA/2
        } else {
          tuningA <- tuningA/1.5
        }
      }
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
    draws[i,1] <- draws[i-1,1]
    draws[i,3] <- draws[i-1,3]
    draws[i,6] <- draws[i-1,6]
    astar <- rnorm(1, draws[i-1,1], tuningA)
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
      lnew <- logPost(c( draws[i,1:5], kstar, draws[i,7:ncol(draws)]), data)
      lold <- logPost(c( draws[i,1:5], draws[i-1,6], draws[i,7:ncol(draws)]), data)
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
    for(j in 7:ncol(finalDraws)){
      d[i] <- d[i] - 
        2*(logMPLPsum(data[[j-6]]$Upper[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 1 & data[[j-6]]$trun == F)], data[[j-6]]$Lower[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 1 & data[[j-6]]$trun == F)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,j]) +
             logMPLPsum(data[[j-6]]$Upper[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 2 & data[[j-6]]$trun == F)], data[[j-6]]$Lower[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 2 & data[[j-6]]$trun == F)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,j]*finalDraws[i,4]) +
             logMPLPsum(data[[j-6]]$Upper[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 3 & data[[j-6]]$trun == F)], data[[j-6]]$Lower[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 3 & data[[j-6]]$trun == F)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,j]*finalDraws[i,4]*finalDraws[i,5]) -
             sum(finalDraws[i,j]*(data[[j-6]]$totalMiles[which(data[[j-6]]$Censor == 1 & data[[j-6]]$Phase == 1 & data[[j-6]]$trun == F)]^finalDraws[i,3])) -
             sum(finalDraws[i,j]*finalDraws[i,4]*(data[[j-6]]$totalMiles[which(data[[j-6]]$Censor == 1 & data[[j-6]]$Phase == 2 & data[[j-6]]$trun == F)]^finalDraws[i,3])) -
             sum(finalDraws[i,j]*finalDraws[i,4]*finalDraws[i,5]*(data[[j-6]]$totalMiles[which(data[[j-6]]$Censor == 1 & data[[j-6]]$Phase == 3 & data[[j-6]]$trun == F)]^finalDraws[i,3])) +
             logMPLPInterval(data[[j-6]]$Lower[which(data[[j-6]]$Phase == 1 & data[[j-6]]$trun == T)], data[[j-6]]$Upper[which(data[[j-6]]$Phase == 1 & data[[j-6]]$trun == T)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,j]) +
             logMPLPInterval(data[[j-6]]$Lower[which(data[[j-6]]$Phase == 2 & data[[j-6]]$trun == T)], data[[j-6]]$Upper[which(data[[j-6]]$Phase == 2 & data[[j-6]]$trun == T)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,j]*finalDraws[i,4]) +
             logMPLPInterval(data[[j-6]]$Lower[which(data[[j-6]]$Phase == 3 & data[[j-6]]$trun == T)], data[[j-6]]$Upper[which(data[[j-6]]$Phase == 3 & data[[j-6]]$trun == T)], finalDraws[i,6], finalDraws[i,3], finalDraws[i,j]*finalDraws[i,4]*finalDraws[i,5]))
    }
  }
  
  davg <- mean(d)
  dthetahat <- 0
  for(j in 7:ncol(finalDraws)){
    dthetahat <- dthetahat - 2*(logMPLPsum(data[[j-6]]$Upper[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 1 & data[[j-6]]$trun == F)], data[[j-6]]$Lower[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 1 & data[[j-6]]$trun == F)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,j])) +
                                  logMPLPsum(data[[j-6]]$Upper[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 2 & data[[j-6]]$trun == F)], data[[j-6]]$Lower[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 2 & data[[j-6]]$trun == F)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,j])*mean(finalDraws[,4])) +
                                  logMPLPsum(data[[j-6]]$Upper[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 3 & data[[j-6]]$trun == F)], data[[j-6]]$Lower[which(data[[j-6]]$Censor == 0 & data[[j-6]]$Phase == 3 & data[[j-6]]$trun == F)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,j])*mean(finalDraws[,4])*mean(finalDraws[,5])) -
                                  sum(mean(finalDraws[,j])*(data[[j-6]]$totalMiles[which(data[[j-6]]$Censor == 1 & data[[j-6]]$Phase == 1 & data[[j-6]]$trun == F)]^mean(finalDraws[,3]))) -
                                  sum(mean(finalDraws[,j])*mean(finalDraws[,4])*(data[[j-6]]$totalMiles[which(data[[j-6]]$Censor == 1 & data[[j-6]]$Phase == 2 & data[[j-6]]$trun == F)]^mean(finalDraws[,3]))) -
                                  sum(mean(finalDraws[,j])*mean(finalDraws[,4])*mean(finalDraws[,5])*(data[[j-6]]$totalMiles[which(data[[j-6]]$Censor == 1 & data[[j-6]]$Phase == 3 & data[[j-6]]$trun == F)]^mean(finalDraws[,3]))) +
                                  logMPLPInterval(data[[j-6]]$Lower[which(data[[j-6]]$Phase == 1 & data[[j-6]]$trun == T)], data[[j-6]]$Upper[which(data[[j-6]]$Phase == 1 & data[[j-6]]$trun == T)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,j])) +
                                  logMPLPInterval(data[[j-6]]$Lower[which(data[[j-6]]$Phase == 2 & data[[j-6]]$trun == T)], data[[j-6]]$Upper[which(data[[j-6]]$Phase == 2 & data[[j-6]]$trun == T)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,j])*mean(finalDraws[,4])) +
                                  logMPLPInterval(data[[j-6]]$Lower[which(data[[j-6]]$Phase == 3 & data[[j-6]]$trun == T)], data[[j-6]]$Upper[which(data[[j-6]]$Phase == 3 & data[[j-6]]$trun == T)], mean(finalDraws[,6]), mean(finalDraws[,3]), mean(finalDraws[,j])*mean(finalDraws[,4])*mean(finalDraws[,5])))
  }
  
  pd <- davg - dthetahat
  dic <- davg + pd
  
  return(list(draws = finalDraws,
              acceptanceShape = accs/samples,
              acceptanceAlpha = acca/samples, 
              acceptanceKappa = acck/samples,
              DIC = dic,
              PD = pd,
              Deviance = d))
}