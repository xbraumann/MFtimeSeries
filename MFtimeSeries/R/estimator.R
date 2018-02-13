#'
#' Wrapper function for the estimation of stable AR and positive semi-definite innovation matrix
#'
#' @param dat.mf T*n data matrix of mf-data
#' @param p Order of the AR polynomial
#' @param method The method for the initialization
#' @param k A predefined k
#' @param upper The upper bound for the k-search
#'
#' @return
#'
#' @export

estimate.AR <- function(dat.mf, p, method="IVL", ...){
  args <- list(...)

  stopifnot(p>0, is.matrix(dat.mf))

  if(method == 'XYW'){
    k <- args[['k']]
    a.hat.0 <- XYW(dat.mf, p=p, k=k)$a.hat
  }

  if(method == 'IVL'){
    k <- args[['k']]
    upper <- args[['upper']]
    a.hat.0 <- IVL(dat.mf, p=p, k=k, upper=upper)$a.hat
    }

  a.hat <- A.make.stable(a.hat.0, ...)$a

  sigma.help <- sigma.hat(dat.mf=dat.mf, a=a.hat, ...)


  list(a=a.hat, a.hat.0=a.hat.0, sigma=sigma.help$sigma, sigma.0=sigma.help$sigma.0)

}


## Mixed-frequency estimators

## XYW

XYW <- function( ... ) {



    }



#'
#' Instrumental Variable (IV) estimator for MF-data
#'
#' @param dat.mf T*n data matrix of mf-data
#' @param p Order of the AR polynomial
#' @param k A predefined k
#' @param upper The upper bound for the k-search
#'
#' @return a.hat The estimator of the AR polynomial
#' @return k The estimated lag
#'
#' @export
IVL <- function(dat.mf, p, k=NULL, upper=NULL) {

  n <- ncol(dat.mf)
  detection <- detect_fast_slow(dat.mf)

  nf <- sum(detection$is.fast)
  ns <- n-nf
  N <- detection$N

  stopifnot(n>1, n>nf, p>0, k >= n*p, N>1)

  if(is.null(k)) k <- AIC.IVL(dat.mf=dat.mf, p=p, n=n, nf=nf, N=N, upper=upper)$k

  dat.f <- matrix( dat.mf[,1:nf], ncol=nf)
  dat.s <- matrix( na.omit(dat.mf[,nf + 1:ns]), ncol=ns)

  x.t.hat <- Yminus(dat.mf, p-1)[1:nrow(dat.mf),]
  Y.minus   <-  Yminus (dat.f, k)

  for(i in 1:p){

    sel <- 1:nrow(dat.mf) + i-1
    Y.minus.move <- Y.minus[sel,]
    coef <- matrix(lm(dat.s ~ 0 + Y.minus.move[c(rep(FALSE, N-1), TRUE), ])$coef, ncol=ns)

    x.t.hat[, n*(i-1) + nf + 1:ns] <- Y.minus[1:nrow(dat.mf),] %*% coef
  }

  a.hat <- matrix(t(lm(x.t.hat[-1, ] ~ 0 + x.t.hat[-nrow(x.t.hat),])$coef)[1:n,],nrow=n)


  list(a.hat=a.hat, k=k)

}



AIC.IVL <- function(dat.mf, p, n, nf, N, upper=NULL){

  ns <- n-nf


  dat.f <- matrix( dat.mf[,1:nf], ncol=nf)
  dat.s <- matrix( na.omit(dat.mf[,nf + 1:ns]), ncol=ns)

  if(is.null(upper)) upper <- log(nrow(dat.s))^1.5
  k.min <- ceiling(n*p / nf)

  stopifnot(upper >= k.min)

  x.t.hat <- Yminus(dat.mf, p-1)[1:nrow(dat.mf),]

  AIC <- NULL

  for (k in k.min:upper){
    Y.minus   <-  Yminus (dat.f, k)

    sel <- 1:nrow(dat.mf)
    Y.minus.move <- Y.minus[sel,]
    resids <- matrix(lm(dat.s ~ 0 + Y.minus.move[c(rep(FALSE, N-1), TRUE), ])$residuals, ncol=ns)

    len <- nrow(resids)

    num.par <- ns * nf * k
    Sigma_ML <- log(abs(det(crossprod(resids)/ len)))

    AIC.help <- Sigma_ML+  2 * num.par / len
    AIC <- c(AIC, AIC.help)

  }
  k <- which(min(AIC) == AIC) + k.min - 1

  list(k=k, AIC=AIC)
}
