library(devtools)
load_all(pkg="./MFtimeSeries")

set.seed(2)

n <- 3
p <- 2
model <- modgen(n=n,p=p)
a <- model$a
b <- model$b
sigma <- model$sigma



Y <- datagen(T=250, a=a, b=b, N=2, nf=2)
dat.mf <- Y$mixed

## IVL
res <- estimate.AR(dat.mf, p)
res$a

## XYW
res <- estimate.AR(dat.mf, p, "XYW", k=n*p)
res$a

a
a-res$a

res$sigma
sigma
sigma-res$sigma
