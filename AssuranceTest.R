


# Want to create test plan that controls for both consumer and producer risk
# Consumer risk: given test is passed, min prob expected failures is large
# P( E(fails in CMiles) > failLim | test passed) < alpha
#
# Producer risk: given test failed, min prob expected failures is small
# P( E(fails in PMiles) < failLim | test failed)
# 
# Cmiles < Pmiles
#
# vehNum: vehicle model number
# n: number of vehicles to test 
# allowedFailes: number of failures allowed to pass test


assuranceTest <- function(postD, vehNum = 1, n = 1, allowedFails = 0,
                          prodRisk = .05, conRisk = .1, failLim = 2,
                          Cmiles = 80, Pmiles = 140,
                          mileStart = 500, increment = 25){
  
  miles <- mileStart
  
  running <- T
  while(running){
    
    if(miles %% 500 == 0){
      numeratorP <- 0
      denominatorP <- 0
      for(j in 1:1500){
        rateParm <- 0
        for(i in 1:20){
          rateParm <- rateParm + 1/postD[[i]]$draws[j,4]*postD[[i]]$draws[j,5]*postD[[i]]$draws[j,vehNum + 5]
        }
        rateParm <- 1/rateParm
        amount <- 0
        for(k in 0:allowedFails){
          amount <- amount + 1 - (exp(-rateParm*(miles^postD[[i]]$draws[j,3])) * (rateParm*(miles^postD[[i]]$draws[j,3]))^k)/factorial(k)
        }
        if((Pmiles^postD[[i]]$draws[j,3])*rateParm < failLim){
          numeratorP <- numeratorP + amount
        }
        denominatorP <- denominatorP + amount
      }
      if(numeratorP/denominatorP > prodRisk){
        print("Producer risk criteria cannot be satisfied")
        print(numeratorP/denominatorP)
        return(NULL)
      }
    }
    
    print(miles)
    
    numeratorC <- 0
    denominatorC <- 0
    for(j in 1:1500){
      rateParm <- 0
      for(i in 1:20){
        rateParm <- rateParm + 1/postD[[i]]$draws[j,4]*postD[[i]]$draws[j,5]*postD[[i]]$draws[j,vehNum + 5]
      }
      rateParm <- 1/rateParm
      amount <- 0
      for(k in 0:allowedFails){
        amount <- amount + (exp(-rateParm*(miles^postD[[i]]$draws[j,3])) * (rateParm*(miles^postD[[i]]$draws[j,3]))^k)/factorial(k)
      }
      if((Cmiles^postD[[i]]$draws[j,3])*rateParm > failLim){
        numeratorC <- numeratorC + amount
      }
      denominatorC <- denominatorC + amount
    }
    
    if(miles %% 500 == 0){
      print("Consumer risk alpha")
      print(numeratorC/denominatorC)
    }
    
    if(numeratorC/denominatorC < conRisk){
      numeratorP <- 0
      denominatorP <- 0
      for(j in 1:1500){
        rateParm <- 0
        for(i in 1:20){
          rateParm <- rateParm + 1/postD[[i]]$draws[j,4]*postD[[i]]$draws[j,5]*postD[[i]]$draws[j,vehNum + 5]
        }
        rateParm <- 1/rateParm
        amount <- 0
        for(k in 0:allowedFails){
          amount <- amount + 1 - (exp(-rateParm*(miles^postD[[i]]$draws[j,3])) * (rateParm*(miles^postD[[i]]$draws[j,3]))^k)/factorial(k)
        }
        if((Pmiles^postD[[i]]$draws[j,3])*rateParm < failLim){
          numeratorP <- numeratorP + amount
        }
        denominatorP <- denominatorP + amount
      }
      
      if(numeratorP/denominatorP < prodRisk){
        running <- F
      } else {
        miles <- miles + increment
      }
    } else {
      miles <- miles + increment
    }
  }
  print("Final test miles")
  return(miles)
}

load("RBAOWeiResults.RData")

assuranceTest(BAOWeiResults, vehNum = 5, failLim = 3, allowedFails = 4,
              mileStart = 1000, increment = 100)



fails1 <- c()

vehNum <- 5

for(j in 1:1500){
  rateParm <- 0
  for(i in 1:20){
    rateParm <- rateParm + BAOWeiResults[[i]]$draws[j,4]*BAOWeiResults[[i]]$draws[j,5]*BAOWeiResults[[i]]$draws[j,vehNum + 5]
  }
  fails1[j] <- (140^BAOWeiResults[[i]]$draws[j,3])*rateParm 
}

hist(fails1)


fails2 <- c()

for(j in 1:1500){
  rateParm <- 0
  for(i in 1:20){
    rateParm <- rateParm + BAOWeiResults[[i]]$draws[j,4]*BAOWeiResults[[i]]$draws[j,5]*BAOWeiResults[[i]]$draws[j,vehNum + 5]
  }
  fails2[j] <- (80^BAOWeiResults[[i]]$draws[j,3])*rateParm 
}

hist(fails2)

plot(density(fails2))
lines(density(fails1))









#################################################################################
#
#################################################################################


assuranceTest <- function(postD, vehNum = 1, n = 1, allowedFails = 0,
                          prodRisk = .05, conRisk = .1, failLim = 2,
                          Cmiles = 80, Pmiles = 140,
                          mileStart = 500, increment = 25){
  
  miles <- mileStart
  
  running <- T
  while(running){
    
    if(miles %% 500 == 0){
      numeratorP <- 0
      denominatorP <- 0
      for(j in 1:1500){
        rateParm <- 0
        for(i in 1:20){
          rateParm <- rateParm + 1/postD[[i]]$draws[j,4]*postD[[i]]$draws[j,5]*postD[[i]]$draws[j,vehNum + 5]
        }
        rateParm <- 1/rateParm
        amount <- 0
        for(k in 0:allowedFails){
          amount <- amount + 1 - (exp(-rateParm*(miles^postD[[i]]$draws[j,3])) * (rateParm*(miles^postD[[i]]$draws[j,3]))^k)/factorial(k)
        }
        if((Pmiles^postD[[i]]$draws[j,3])*rateParm < failLim){
          numeratorP <- numeratorP + amount
        }
        denominatorP <- denominatorP + amount
      }
      if(numeratorP/denominatorP > prodRisk){
        print("Producer risk criteria cannot be satisfied")
        print(numeratorP/denominatorP)
        return(NULL)
      }
    }
    
    print(miles)
    
    numeratorC <- 0
    denominatorC <- 0
    for(j in 1:1500){
      rateParm <- 0
      for(i in 1:20){
        rateParm <- rateParm + 1/postD[[i]]$draws[j,4]*postD[[i]]$draws[j,5]*postD[[i]]$draws[j,vehNum + 5]
      }
      rateParm <- 1/rateParm
      amount <- 0
      for(k in 0:allowedFails){
        amount <- amount + (exp(-rateParm*(miles^postD[[i]]$draws[j,3])) * (rateParm*(miles^postD[[i]]$draws[j,3]))^k)/factorial(k)
      }
      if((Cmiles^postD[[i]]$draws[j,3])*rateParm > failLim){
        numeratorC <- numeratorC + amount
      }
      denominatorC <- denominatorC + amount
    }
    
    if(miles %% 500 == 0){
      print("Consumer risk alpha")
      print(numeratorC/denominatorC)
    }
    
    if(numeratorC/denominatorC < conRisk){
      numeratorP <- 0
      denominatorP <- 0
      for(j in 1:1500){
        rateParm <- 0
        for(i in 1:20){
          rateParm <- rateParm + 1/postD[[i]]$draws[j,4]*postD[[i]]$draws[j,5]*postD[[i]]$draws[j,vehNum + 5]
        }
        rateParm <- 1/rateParm
        amount <- 0
        for(k in 0:allowedFails){
          amount <- amount + 1 - (exp(-rateParm*(miles^postD[[i]]$draws[j,3])) * (rateParm*(miles^postD[[i]]$draws[j,3]))^k)/factorial(k)
        }
        if((Pmiles^postD[[i]]$draws[j,3])*rateParm < failLim){
          numeratorP <- numeratorP + amount
        }
        denominatorP <- denominatorP + amount
      }
      
      if(numeratorP/denominatorP < prodRisk){
        running <- F
      } else {
        miles <- miles + increment
      }
    } else {
      miles <- miles + increment
    }
  }
  print("Final test miles")
  return(miles)
}

load("RWeibullResults.RData")

assuranceTest(WeibullResults, vehNum = 1, failLim = 3, allowedFails = 2,
              mileStart = 1000, increment = 100)



fails1 <- c()

for(j in 1:1500){
  rateParm <- 0
  for(i in 1:20){
    rateParm <- rateParm + WeibullResults[[i]]$draws[j,4]*WeibullResults[[i]]$draws[j,5]*WeibullResults[[i]]$draws[j,6]
  }
  fails1[j] <- (140^WeibullResults[[i]]$draws[j,3])*rateParm 
}

hist(fails1)


fails2 <- c()

for(j in 1:1500){
  rateParm <- 0
  for(i in 1:20){
    rateParm <- rateParm + WeibullResults[[i]]$draws[j,4]*WeibullResults[[i]]$draws[j,5]*WeibullResults[[i]]$draws[j,6]
  }
  fails2[j] <- (80^WeibullResults[[i]]$draws[j,3])*rateParm 
}

hist(fails2)

plot(density(fails2))
lines(density(fails1))





#################################################################################
#
#################################################################################



assuranceTest <- function(postD, vehNum = 1, n = 1, allowedFails = 0,
                          prodRisk = .05, conRisk = .1, failLim = 2,
                          Cmiles = 80, Pmiles = 140,
                          mileStart = 500, increment = 25){
  
  miles <- mileStart
  
  running <- T
  while(running){
    
    if(miles %% 500 == 0){
      numeratorP <- 0
      denominatorP <- 0
      for(j in 1:1500){
        rateParm <- 0
        for(i in 1:20){
          rateParm <- rateParm + postD[[i]]$draws[j,4]*postD[[i]]$draws[j,5]*postD[[i]]$draws[j,vehNum + 5]
        }
        amount <- 0
        for(k in 0:allowedFails){
          amount <- amount + 1 - (exp(-rateParm*(miles^postD[[i]]$draws[j,3])) * (rateParm*(miles^postD[[i]]$draws[j,3]))^k)/factorial(k)
        }
        if((Pmiles^postD[[i]]$draws[j,3])*rateParm < failLim){
          numeratorP <- numeratorP + amount
        }
        denominatorP <- denominatorP + amount
      }
      if(numeratorP/denominatorP > prodRisk){
        print("Producer risk criteria cannot be satisfied")
        print(numeratorP/denominatorP)
        return(NULL)
      }
    }
    
    print(miles)
    
    numeratorC <- 0
    denominatorC <- 0
    for(j in 1:1500){
      rateParm <- 0
      for(i in 1:20){
        rateParm <- rateParm + postD[[i]]$draws[j,4]*postD[[i]]$draws[j,5]*postD[[i]]$draws[j,vehNum + 5]
      }
      amount <- 0
      for(k in 0:allowedFails){
        amount <- amount + (exp(-rateParm*(miles^postD[[i]]$draws[j,3])) * (rateParm*(miles^postD[[i]]$draws[j,3]))^k)/factorial(k)
      }
      if((Cmiles^postD[[i]]$draws[j,3])*rateParm > failLim){
        numeratorC <- numeratorC + amount
      }
      denominatorC <- denominatorC + amount
    }
    
    if(miles %% 500 == 0){
      print("Consumer risk alpha")
      print(numeratorC/denominatorC)
    }
    
    if(numeratorC/denominatorC < conRisk){
      numeratorP <- 0
      denominatorP <- 0
      for(j in 1:1500){
        rateParm <- 0
        for(i in 1:20){
          rateParm <- rateParm + 1/postD[[i]]$draws[j,4]*postD[[i]]$draws[j,5]*postD[[i]]$draws[j,vehNum + 5]
        }
        rateParm <- 1/rateParm
        amount <- 0
        for(k in 0:allowedFails){
          amount <- amount + 1 - (exp(-rateParm*(miles^postD[[i]]$draws[j,3])) * (rateParm*(miles^postD[[i]]$draws[j,3]))^k)/factorial(k)
        }
        if((Pmiles^postD[[i]]$draws[j,3])*rateParm < failLim){
          numeratorP <- numeratorP + amount
        }
        denominatorP <- denominatorP + amount
      }
      
      if(numeratorP/denominatorP < prodRisk){
        running <- F
      } else {
        miles <- miles + increment
      }
    } else {
      miles <- miles + increment
    }
  }
  print("Final test miles")
  return(miles)
}

load("RBAOWeiResults.RData")

assuranceTest(BAOWeiResults, vehNum = 5, failLim = 1, allowedFails = 4,
              mileStart = 100, increment = 100)



fails1 <- c()

vehNum <- 5

for(j in 1:1500){
  rateParm <- 0
  for(i in 1:20){
    rateParm <- rateParm + BAOWeiResults[[i]]$draws[j,4]*BAOWeiResults[[i]]$draws[j,5]*BAOWeiResults[[i]]$draws[j,vehNum + 5]
  }
  fails1[j] <- (140^BAOWeiResults[[i]]$draws[j,3])*rateParm 
}

hist(fails1)


fails2 <- c()

for(j in 1:1500){
  rateParm <- 0
  for(i in 1:20){
    rateParm <- rateParm + BAOWeiResults[[i]]$draws[j,4]*BAOWeiResults[[i]]$draws[j,5]*BAOWeiResults[[i]]$draws[j,vehNum + 5]
  }
  fails2[j] <- (80^BAOWeiResults[[i]]$draws[j,3])*rateParm 
}

hist(fails2)

plot(density(fails2))
lines(density(fails1))