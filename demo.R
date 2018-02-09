library(devtools)
load_all(pkg="../MFtimeSeries")

set.seed(1)

model <- modgen(n=2,p=2)
a <- model$a
sigma <- model$sigma



Y <- datagen(T=25, a=a, b=t(chol(sigma)), N=2, nf=1)

str(Y)

head(Y$full)

head(Y$mf2)

acf(Y$full)
pacf(Y$full)
