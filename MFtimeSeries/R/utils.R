

A.companionform <- function(a) {
  stopifnot(is.matrix(a))

  n <- nrow(a)
  np <- ncol(a)

  A <- rbind(a, cbind(diag(np-n),matrix(0, nrow=(np-n), ncol=n)))
  A
}

A.stable <- function(A) {
  if(nrow(A)!=ncol(A)) A <- A.companionform(A)
  return(max(abs(eigen(A, only.values=TRUE)$values)) < 1)
}

calc.acf <- function(data, lags=0){
  stats::acf(data, na.action=na.pass, type="covariance", lag.max=lags, demean=TRUE, plot=FALSE)$acf
}

detect_fast_slow <- function (x) {

  tmp <- apply(x, 2, function(x) {any(is.na(x)) })

  grid <- which(apply(x[1:100,], 1, function(x) {any(is.na(x)) }))
  N <- grid[2]-grid[1]

  list(is.fast = !tmp,
       is.slow = tmp,
       N = N)
}


Yminus <- function (data, k) {

  d <- ncol(data)
  T <- nrow(data)

  Y <- NULL
  for (i in 0:k) {
      # unavailable values are set to 0
      ith.col <- rbind(matrix(0, ncol=d,  nrow=i),
                       data,
                       matrix(0, ncol=d,  nrow=k-i))

      Y <- cbind(Y, ith.col)
  }
  Y
}

