library(devtools)
load_all(pkg="../MFtimeSeries")

n <- 2
p <- 2

set.seed(1)

model <- modgen(n,p)
a <- model$a
sigma <- model$sigma
