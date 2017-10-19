using Distributions
using DataFrames
using CSV

α = .4
β = .1

lambdas = rand(Gamma(α, β), 8)

data = zeros(500, 8)

for i in 1:8
  data[:,i] = rand(Exponential(lambdas[i]), 500)
end

CSV.write("simExpData2.csv", convert(DataFrame, data))
