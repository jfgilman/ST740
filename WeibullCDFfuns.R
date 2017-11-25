

# CDF of weibull
wF <- function(x, scale, shape){
  return(1 - exp(-scale*(x)^shape))
}

# inverse CDF of weibull for generating data
wF.inv <- function(u,lam,beta){
  return((-log(1-u)/lam)^(1/beta))
}

plpF <- function(y, start, scale, shape){
  return(exp(scale*(start)^shape) * (1 - exp(-scale*(y)^shape)))
}

plpF.inv <- function(u, start, scale, shape){
  return((-log(exp(-scale*(start)^shape)-(u/exp(scale*(start)^shape)))/scale)^(1/shape))
}