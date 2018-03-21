context('Test if estimation procedures give stable estimates')


#library(devtools)
#load_all(pkg="./MFtimeSeries", recompile=TRUE)

## generate a model

m <- modgen(n=2, p=1)

abs(eigen(m$a)$values)



