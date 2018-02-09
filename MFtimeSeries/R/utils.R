

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

  grid <- which(apply(x, 1, function(x) {any(is.na(x)) }))
  N <- grid[2]-grid[1]

  list(is.fast = !tmp,
       is.slow = tmp,
       N = N)
}


#open
Yplus <- function (data, k) {

  d <- dim(data)[2]
  T <- dim(data)[1]
  Y <- NULL
  if(k>0){
    for (i in 1:k) {

      zeros <- matrix(0, nrow=d, ncol=i)
      sel <- i:T
      ith.row <- cbind(t(data)[, sel, drop=FALSE], zeros)

      Y <- rbind(Y, ith.row)
    }
  }
  Y
}

#open
Yminus <- function (data, k, offset=0) {

  d <- dim(data)[2]
  T <- dim(data)[1]
  Y <- NULL
  if(k>0){
    for (i in 1:k) {
      # unavailable values are set to 0
      zerocols <- (offset-i+1<=0)*abs(offset-i)
      zeros <- matrix(0, nrow=d, ncol=zerocols)
      sel <- (offset-i+1):(T-i+1)
      sel <- sel[sel>0]
      ith.row <- cbind(zeros, t(data)[, sel, drop=FALSE])

      Y <- rbind(Y, ith.row)
    }
  }
  Y
}

