library(devtools)
load_all(pkg="../MFtimeSeries")

n <- 2
p <- 2

set.seed(1)

model <- modgen(n,p)
a <- model$a
sigma <- model$sigma



Y <- datagen(25, a, t(chol(sigma)), 2, 1)

str(Y)

head(Y$full)

head(Y$mf2)

acf(Y$full)
pacf(Y$full)
