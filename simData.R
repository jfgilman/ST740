
# Create Simulated data

alpha <- .4
beta <- .1

lambdas <- rgamma(8, alpha, beta)

simData <- c()

for(i in 1:length(lambdas)){
  simData <- cbind(simData, rexp(500, lambdas[i]))
}

write.csv(simData, "simExpData.csv")
