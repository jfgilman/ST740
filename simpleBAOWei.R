library(Rcpp)

cppFunction('double logWeibullSumCPP(NumericVector x, double shape, double scale) {
            
            int nx = x.size();
            double d = 0.0;
            
            for(int i = 0; i < nx; i++) {
            d = d + log(shape) + log(scale) + (shape - 1)*log(x[i]) - scale*pow(x[i],shape);
            }
            
            return d;
            }')

simpBAOWeiMCMC <- function(data, samples = 5000, shapePriorA = .001,
                        shapePriorB = .001, lamPriorA = .001, lamPriorB = .001,
                        theta1PriorA = .001, theta1PriorB = .001, theta2PriorA = .001,
                        theta2PriorB = .001, theta1Start = 1, theta2Start = 1, betaStart = 1,
                        betaPriorA = .001, betaPriorB = .001, tuning = 1, 
                        burnin = 1000, thin = 10, tuningS = 1){
  
  # matrix for keeping MCMC draws for each parameter
  lam_draws <- rep(0, samples)
  lam_draws[1] <- 1
  
  beta_draws <- rep(0, samples)
  beta_draws[1] <- betaStart
  
  theta1_draws <- rep(0, samples)
  theta1_draws[1] <- theta1Start
  
  theta2_draws <- rep(0, samples)
  theta2_draws[1] <- theta2Start
  
  accs <- 0
  
  # log posterior function for MCMC draws
  logPost <- function(parm, d) {
    
    lp <- 0
    for(k in 1:length(d)){
      lp <- lp + 
        sum(logWeibullSumCPP(d[[k]]$totalMiles[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 1 & d[[k]]$trun == F)], parm[4], parm[1])) +
        sum(logWeibullSumCPP(d[[k]]$totalMiles[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 2 & d[[k]]$trun == F)], parm[4], parm[1]*parm[2])) +
        sum(logWeibullSumCPP(d[[k]]$totalMiles[which(d[[k]]$Censor == 0 & d[[k]]$Phase == 3 & d[[k]]$trun == F)], parm[4], parm[1]*parm[2]*parm[3])) -
        sum(parm[1]*(d[[k]]$totalMiles[which(d[[k]]$Censor == 1 & d[[k]]$Phase == 1 & d[[k]]$trun == F)]^parm[4])) -
        sum(parm[1]*parm[2]*(d[[k]]$totalMiles[which(d[[k]]$Censor == 1 & d[[k]]$Phase == 2 & d[[k]]$trun == F)]^parm[4])) -
        sum(parm[1]*parm[2]*parm[3]*(d[[k]]$totalMiles[which(d[[k]]$Censor == 1 & d[[k]]$Phase == 3 & d[[k]]$trun == F)]^parm[4])) +
        sum(parm[1]*(d[[k]]$Upper[which(d[[k]]$Phase == 1 & d[[k]]$trun == T)]^parm[4])) +
        sum(parm[1]*parm[2]*(d[[k]]$Upper[which(d[[k]]$Phase == 2 & d[[k]]$trun == T)]^parm[4])) + 
        sum(parm[1]*parm[2]*parm[3]*(d[[k]]$Upper[which(d[[k]]$Phase == 3 & d[[k]]$trun == T)]^parm[4])) 
    }
    
    lp <- lp + 
      dgamma(parm[1], lamPriorA, lamPriorB, log=T) + 
      dgamma(parm[2], theta1PriorA, theta1PriorB, log=T) +
      dgamma(parm[3], theta2PriorA, theta2PriorB, log=T) +
      dgamma(parm[4], betaPriorA, betaPriorB, log=T)
    
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
  
  # MCMC draws
  for (i in 2:samples) {
    
    datSum <- 0
    for(j in 1:length(data)){
      datSum <- datSum + sum(data[[j]]$totalMiles[which(data[[j]]$Phase == 1)],
                             (theta1_draws[i-1] * data[[j]]$totalMiles[which(data[[j]]$Phase == 2)]),
                             (theta1_draws[i-1] * theta2_draws[i-1] * data[[j]]$totalMiles[which(data[[j]]$Phase == 3)]))
    }
    
    lam_draws[i] <- rgamma(1, noCenCount + lamPriorA, datSum + lamPriorB)
    
    # Phase sums
    phase2Sum <- 0
    phase3Sum4 <- 0
    phase3Sum5 <- 0
    for(k in 1:length(data)){
      phase2Sum <- phase2Sum + lam_draws[i] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 2)])
      
      # for the first theta
      phase3Sum4 <- phase3Sum4 + theta2_draws[i-1] * lam_draws[i] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 3)])
      
      # for second theta
      phase3Sum5 <- phase3Sum5 + theta1_draws[i-1] * lam_draws[i] * sum(data[[k]]$totalMiles[which(data[[k]]$Phase == 3)])
    }
    
    theta1_draws[i] <- rgamma(1, phase2Count + phase3Count + theta1PriorA, phase2Sum + phase3Sum4 + theta1PriorB)
    theta2_draws[i] <- rgamma(1, phase3Count + theta2PriorA, phase3Sum5 + theta2PriorB)
    
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
    
    beta_draws[i] <- beta_draws[i-1]
    # Sample from shape
    sstar <- rnorm(1, beta_draws[i-1], tuningS)
    if (sstar > 0) {
      lnew <- logPost(c(lam_draws[i], theta1_draws[i], theta2_draws[i], sstar), data)
      lold <- logPost(c(lam_draws[i], theta1_draws[i], theta2_draws[i], beta_draws[i-1]), data)
      if(is.finite(lnew - lold)){
        if (lnew - lold > log(runif(1))) {
          beta_draws[i] <- sstar
          accs <- accs + 1
        }
      }
    }
  }
  
  lam_finalDraws <- lam_draws[seq(from = burnin + 1, to = samples, by = thin)]
  theta1_finalDraws <- theta1_draws[seq(from = burnin + 1, to = samples, by = thin)]
  theta2_finalDraws <- theta2_draws[seq(from = burnin + 1, to = samples, by = thin)]
  beta_finalDraws <- beta_draws[seq(from = burnin + 1, to = samples, by = thin)]
  
  # DIC
  d <- rep(0, length(lam_finalDraws))
  for(i in 1:length(lam_finalDraws)){
    for(j in 1:length(data)){
      d[i] <- d[i] - 
        2*(sum(logWeibullSumCPP(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1 & data[[j]]$trun == F)], beta_finalDraws[i], lam_finalDraws[i])) +
             sum(logWeibullSumCPP(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2 & data[[j]]$trun == F)], beta_finalDraws[i], lam_finalDraws[i]*theta1_finalDraws[i])) +
             sum(logWeibullSumCPP(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3 & data[[j]]$trun == F)], beta_finalDraws[i], lam_finalDraws[i]*theta1_finalDraws[i]*theta2_finalDraws[i])) -
             sum(lam_finalDraws[i]*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 1 & data[[j]]$trun == F)]^beta_finalDraws[i])) -
             sum(lam_finalDraws[i]*theta1_finalDraws[i]*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 2 & data[[j]]$trun == F)]^beta_finalDraws[i])) -
             sum(lam_finalDraws[i]*theta1_finalDraws[i]*theta2_finalDraws[i]*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 3 & data[[j]]$trun == F)]^beta_finalDraws[i])) +
             sum(lam_finalDraws[i]*(data[[j]]$Upper[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)]^beta_finalDraws[i])) +
             sum(lam_finalDraws[i]*theta1_finalDraws[i]*(data[[j]]$Upper[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)]^beta_finalDraws[i])) + 
             sum(lam_finalDraws[i]*theta1_finalDraws[i]*theta2_finalDraws[i]*(data[[j]]$Upper[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)]^beta_finalDraws[i])))
    }
  }
  
  davg <- mean(d)
  dthetahat <- 0
  for(j in 1:length(data)){
    dthetahat <- dthetahat - 2*(sum(logWeibullSumCPP(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 1 & data[[j]]$trun == F)], mean(beta_finalDraws), mean(lam_finalDraws))) +
                                  sum(logWeibullSumCPP(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 2 & data[[j]]$trun == F)], mean(beta_finalDraws), mean(lam_finalDraws)*mean(theta1_finalDraws))) +
                                  sum(logWeibullSumCPP(data[[j]]$totalMiles[which(data[[j]]$Censor == 0 & data[[j]]$Phase == 3 & data[[j]]$trun == F)], mean(beta_finalDraws), mean(lam_finalDraws)*mean(theta1_finalDraws)*mean(theta2_finalDraws))) -
                                  sum(mean(lam_finalDraws)*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 1 & data[[j]]$trun == F)]^mean(beta_finalDraws))) -
                                  sum(mean(lam_finalDraws)*mean(theta1_finalDraws)*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 2 & data[[j]]$trun == F)]^mean(beta_finalDraws))) -
                                  sum(mean(lam_finalDraws)*mean(theta1_finalDraws)*mean(theta2_finalDraws)*(data[[j]]$totalMiles[which(data[[j]]$Censor == 1 & data[[j]]$Phase == 3 & data[[j]]$trun == F)]^mean(beta_finalDraws))) +
                                  sum(mean(lam_finalDraws)*(data[[j]]$Upper[which(data[[j]]$Phase == 1 & data[[j]]$trun == T)]^mean(beta_finalDraws))) +
                                  sum(mean(lam_finalDraws)*mean(theta1_finalDraws)*(data[[j]]$Upper[which(data[[j]]$Phase == 2 & data[[j]]$trun == T)]^mean(beta_finalDraws))) + 
                                  sum(mean(lam_finalDraws)*mean(theta1_finalDraws)*mean(theta2_finalDraws)*(data[[j]]$Upper[which(data[[j]]$Phase == 3 & data[[j]]$trun == T)]^mean(beta_finalDraws))))
  }
  
  pd <- davg - dthetahat
  dic <- davg + pd
  
  return(list(lam_draws = lam_finalDraws,
              theta1_draws = theta1_finalDraws,
              theta2_draws = theta2_finalDraws,
              shape_draws = beta_finalDraws,
              DIC = dic,
              PD = pd,
              Deviance = d))
}