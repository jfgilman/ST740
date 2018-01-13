

m0 = 62

dpois(0, lambda = (m0*(3/(80^1))))

1 - dpois(0, lambda = (m0*(3/(140^1))))


m1 = 104

dpois(0, lambda = (m1*(3/(80^1)))) + dpois(1, lambda = (m1*(3/(80^1))))

1 - dpois(0, lambda = (m1*(3/(140^1)))) - dpois(1, lambda = (m1*(3/(140^1))))


m2 = 142

dpois(0, lambda = (m2*(3/(80^1)))) + dpois(1, lambda = (m2*(3/(80^1)))) + dpois(2, lambda = (m2*(3/(80^1))))

1 - dpois(0, lambda = (m2*(3/(140^1)))) - dpois(1, lambda = (m2*(3/(140^1)))) - dpois(2, lambda = (m2*(3/(140^1))))


m3 = 179

dpois(0, lambda = (m3*(3/(80^1)))) + dpois(1, lambda = (m3*(3/(80^1)))) + dpois(2, lambda = (m3*(3/(80^1)))) + dpois(3, lambda = (m3*(3/(80^1))))

1 - dpois(0, lambda = (m3*(3/(140^1)))) - dpois(1, lambda = (m3*(3/(140^1)))) - dpois(2, lambda = (m3*(3/(140^1)))) - dpois(3, lambda = (m3*(3/(140^1))))


freqExpTest <- function(aFails = 0, cRate = 3/(80^.3), pRate = 3/(140^.3), startM = 50){
  prob1 <- 1
  miles <- startM - 5
  while(prob1 > 0.05){
    prob1 <- 0
    miles <- miles + 5
    for(i in 0:aFails){
      prob1 = prob1 + dpois(i, lambda = (miles*cRate))
    }
  }
  
  prob2 <- 1
  for(i in 0:aFails){
    prob2 = prob2 - dpois(i, lambda = (miles*pRate))
  }

  return(c(miles, prob1, prob2))
}

for(i in 200:300){
  print(freqExpTest(i))
}
