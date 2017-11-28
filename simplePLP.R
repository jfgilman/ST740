# library(Rcpp)
# 
# cppFunction('double logPLPSum(NumericVector x, double shape, double scale) {
# 
#             int nx = x.size();
#             double d = 0.0;
# 
#             for(int i = 0; i < nx; i++) {
#             d = d + log(shape) + log(scale) + (shape - 1)*log(x[i]);
#             }
# 
#             return d;
#             }')
# 
# cppFunction('double logPLPInterval(NumericVector L, NumericVector U, double shape, double scale) {
# 
#             int nx = L.size();
#             double d = 0.0;
# 
#             for(int i = 0; i < nx; i++) {
#             d = d + log(scale) + log(pow(U[i],shape) - pow(L[i],shape)) - scale*(pow(U[i],shape) - pow(L[i],shape));
#             }
# 
#             return d;
#             }')

logPLPSum <- function(x, shape, scale){
  d <- 0
  if(length(x) > 0){
    for(i in 1:length(x)){
      d <- d + log(shape) + log(scale) + (shape - 1)*log(x[i])
    }
  }
  return(d)
}

logPLPInterval <- function(L, U, shape, scale){
  d <- 0
  if(length(L) > 0){
    for(i in 1:length(L)){
      d <- sum(log(scale) + log(U[i]^shape - L[i]^shape) - scale*(U[i]^shape - L[i]^shape))
    }
  }
  return(d)
}



simpPLP <- function(data, samples = 40000, shapePriorA = .001,
                    shapePriorB = .001, hyperT1A = 3, hyperT1B = 3, hyperT2A = 3, hyperT2B = 3,
                    theta1Start = 1, theta2Start = 1, shapeStart = 1, tuningA = 1,
                    tuningS = 1, burnin = 20000, thin = 10, lamPriorA = .001, lamPriorB = .001){
  
  
  # matrix for keeping MCMC draws for each parameter
  draws <- matrix(0, nrow = samples, ncol = 6)

  draws[1,3] <- shapeStart
  draws[1,4] <- theta1Start
  draws[1,5] <- theta2Start
  
  # counter for acceptance rate
  accs <- 0
  
  # log posterior function for MCMC draws
  logPost <- function(parm, d) {
    
    lp <- 0
    
    for(k in 1:length(d)){
      
      lp <- lp + 
        logPLPSum(d[[k]]$totalMiles[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 1 & d[[k]]$trun == F)], parm[3], parm[6]) +
        logPLPSum(d[[k]]$totalMiles[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 2 & d[[k]]$trun == F)], parm[3], parm[6]*parm[4]) +
        logPLPSum(d[[k]]$totalMiles[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 3 & d[[k]]$trun == F)], parm[3], parm[6]*parm[4]*parm[5]) -
        sum(parm[6]*(d[[k]]$totalMiles[which(d[[k]]$Censor == 1 & d[[k]]$Phase == 1 & d[[k]]$trun == F)]^parm[3])) -
        sum(parm[6]*parm[4]*(d[[k]]$totalMiles[which(d[[k]]$Censor == 1 & d[[k]]$Phase == 2 & d[[k]]$trun == F)]^parm[3])) -
        sum(parm[6]*parm[4]*parm[5]*(d[[k]]$totalMiles[which(d[[k]]$Censor == 1 & d[[k]]$Phase == 3 & d[[k]]$trun == F)]^parm[3])) +
        logPLPInterval(d[[k]]$Lower[which(d[[k]]$Phase == 1 & d[[k]]$trun == T)], d[[k]]$Upper[which(d[[k]]$Phase == 1 & d[[k]]$trun == T)], parm[3], parm[6]) +
        logPLPInterval(d[[k]]$Lower[which(d[[k]]$Phase == 2 & d[[k]]$trun == T)], d[[k]]$Upper[which(d[[k]]$Phase == 2 & d[[k]]$trun == T)], parm[3], parm[6]*parm[4]) +
        logPLPInterval(d[[k]]$Lower[which(d[[k]]$Phase == 3 & d[[k]]$trun == T)], d[[k]]$Upper[which(d[[k]]$Phase == 3 & d[[k]]$trun == T)], parm[3], parm[6]*parm[4]*parm[5])
    }
    lp <- lp + 
      dgamma(parm[3], shapePriorA, shapePriorB, log=T) +
      dgamma(parm[4],hyperT1A, hyperT1B, log=T) +
      dgamma(parm[5],hyperT2A, hyperT2B, log=T) +
      dgamma(parm[6], lamPriorA, lamPriorB, log=T)
    
    return(lp)
  }
  
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
  
  for (i in 2:samples) {
    
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
    
    draws[i,6] <- rgamma(1, noCenCount + lamPriorA, datSum + lamPriorB)
    
    # Phase sums
    phase2Sum <- 0
    phase3Sum4 <- 0
    phase3Sum5 <- 0
    for(k in 1:length(data)){
      phase2Sum <- phase2Sum + draws[i,6] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 2 & data[[k]]$Censor == 1 & data[[k]]$trun == F)]^draws[i-1,3],
                                                  data[[k]]$Upper[which(data[[k]]$Phase == 2 & data[[k]]$trun == T)]^draws[i-1,3] -
                                                    data[[k]]$Lower[which(data[[k]]$Phase == 2 & data[[k]]$trun == T)]^draws[i-1,3])
      
      # for the first theta
      phase3Sum4 <- phase3Sum4 + draws[i,6] * draws[i-1,5] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 3 & data[[k]]$Censor == 1 & data[[k]]$trun == F)]^draws[i-1,3],
                                                                   data[[k]]$Upper[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3] -
                                                                     data[[k]]$Lower[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3])
      
      # for second theta
      phase3Sum5 <- phase3Sum5 + draws[i,6] * draws[i-1,4] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 3 & data[[k]]$Censor == 1 & data[[k]]$trun == F)]^draws[i-1,3],
                                                                   data[[k]]$Upper[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3] -
                                                                     data[[k]]$Lower[which(data[[k]]$Phase == 3 & data[[k]]$trun == T)]^draws[i-1,3])
    }
    
    draws[i,4] <- rgamma(1, phase2Count + phase3Count + hyperT1A, phase2Sum + phase3Sum4 + hyperT1B)
    draws[i,5] <- rgamma(1, phase3Count + hyperT2A, phase3Sum5 + hyperT2B)
    
    
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
    }
    
    draws[i,3] <- draws[i-1,3]
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
  }
  
  finalDraws <- draws[seq(from = burnin + 1, to = samples, by = thin),]
  
  # DIC
  d <- rep(0, nrow(finalDraws))
  for(i in 1:nrow(finalDraws)){
    for(j in 1:length(data)){
      d[i] <- d[i] - 
        2*(logPLPSum(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1 & data[[j]]$trun == F)], finalDraws[i,3], finalDraws[i,6]) +
             logPLPSum(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2 & data[[j]]$trun == F)], finalDraws[i,3], finalDraws[i,6]*finalDraws[i,4]) +
             logPLPSum(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3 & data[[j]]$trun == F)], finalDraws[i,3], finalDraws[i,6]*finalDraws[i,4]*finalDraws[i,5]) -
             sum(finalDraws[i,6]*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 1 & data[[j]]$trun == F)]^finalDraws[i,3])) -
             sum(finalDraws[i,6]*finalDraws[i,4]*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 2 & data[[j]]$trun == F)]^finalDraws[i,3])) -
             sum(finalDraws[i,6]*finalDraws[i,4]*finalDraws[i,5]*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 3 & data[[j]]$trun == F)]^finalDraws[i,3])) +
             logPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)], finalDraws[i,3], finalDraws[i,6]) +
             logPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)], finalDraws[i,3], finalDraws[i,6]*finalDraws[i,4]) +
             logPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)], finalDraws[i,3], finalDraws[i,6]*finalDraws[i,4]*finalDraws[i,5]))
    }
  }
  
  davg <- mean(d)
  dthetahat <- 0
  for(j in 1:length(data)){
    dthetahat <- dthetahat - 2*(logPLPSum(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1 & data[[j]]$trun == F)], mean(finalDraws[,3]), mean(finalDraws[,6])) +
                                  logPLPSum(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2 & data[[j]]$trun == F)], mean(finalDraws[,3]), mean(finalDraws[,6])*mean(finalDraws[,4])) +
                                  logPLPSum(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3 & data[[j]]$trun == F)], mean(finalDraws[,3]), mean(finalDraws[,6])*mean(finalDraws[,4])*mean(finalDraws[,5])) -
                                  sum(mean(finalDraws[,6])*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 1 & data[[j]]$trun == F)]^mean(finalDraws[,3]))) -
                                  sum(mean(finalDraws[,6])*mean(finalDraws[,4])*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 2 & data[[j]]$trun == F)]^mean(finalDraws[,3]))) -
                                  sum(mean(finalDraws[,6])*mean(finalDraws[,4])*mean(finalDraws[,5])*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 3 & data[[j]]$trun == F)]^mean(finalDraws[,3]))) +
                                  logPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)], mean(finalDraws[,3]), mean(finalDraws[,6])) +
                                  logPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)], mean(finalDraws[,3]), mean(finalDraws[,6])*mean(finalDraws[,4])) +
                                  logPLPInterval(data[[j]]$Lower[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)], data[[j]]$Upper[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)], mean(finalDraws[,3]), mean(finalDraws[,6])*mean(finalDraws[,4])*mean(finalDraws[,5])))
  }
  
  pd <- davg - dthetahat
  dic <- davg + pd
  
  return(list(draws = finalDraws,
              acceptanceShape = accs/samples,
              DIC = dic,
              PD = pd,
              Deviance = d))
}