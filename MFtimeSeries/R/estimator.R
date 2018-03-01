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

## XYW - Extended Yule Walker

XYW <- function( ... ) {



}



#'
#' Instrumental Variable (IVL) estimator for MF-data
#'
#' The method is described in L. Koelbl, M. Deistler (2018): A New Approach for Estimating VAR Systems in the Mixed-Frequency Case (Statistical Papers)
#' The main idea of this new procedure is to project the slow components on the present and past fast ones in order to create instrumental variables.
#'
#' @param dat.mf The mixed-frequency data (T*n matrix), where the first nf components are the fast ones and the last ns components are the slow.
#' @param p Order of the AR polynomial
#' @param k A predefined k. If it is NULL, it will be estimated
#' @param upper The upper bound for the k estimation. If it is NULL, then log(n_s)^1.5 will be used where n_s is the amount of observed slow components.
#'
#' @return a.hat The estimator of the AR polynomial
#' @return k The estimated lag for the
#'
#' @export
IVL <- function(dat.mf, p, k=NULL, upper=NULL) {

  n <- ncol(dat.mf)
  detection <- detect_fast_slow(dat.mf)

  nf <- sum(detection$is.fast)
  ns <- n-nf
  N <- detection$N
  grid <- !is.na(dat.mf[,nf+1])

  stopifnot(n>1, n>nf, p>0, k>=n*p, N>1)

  if(is.null(k)) k <- AIC.IVL(dat.mf=dat.mf, p=p, n=n, nf=nf, N=N, grid=grid, upper=upper)$k

  dat.f <- dat.mf[,1:nf, drop=FALSE]
  dat.s <- na.omit(dat.mf[,nf + 1:ns, drop=FALSE])

  x.t.hat <- Yminus(dat.mf, p-1)[1:nrow(dat.mf),]
  Y.minus <- Yminus (dat.f, k)

  for(i in 1:p){

    sel <- 1:nrow(dat.mf) + i-1
    Y.minus.move <- Y.minus[sel,]
    coef <- matrix(lm(dat.s ~ 0 + Y.minus.move[grid, ])$coef, ncol=ns)

    x.t.hat[, n*(i-1) + nf + 1:ns] <- Y.minus[1:nrow(dat.mf),] %*% coef
  }

  if(p>1){
    x.t.hat <- x.t.hat[-(1:(p-1)),]
  }


  A <- t(lm(  x.t.hat[-1, ] ~ 0 + x.t.hat[-nrow(x.t.hat),]  )$coef)
  a.hat <- matrix(A[1:n,],nrow=n)

  list(a.hat=a.hat, k=k)

}


#' @describeIn IVL estimation of the k parameter
AIC.IVL <- function(dat.mf, p, n, nf, N, grid, upper=NULL){

  ns <- n-nf

  dat.f <- dat.mf[,1:nf, drop=FALSE]
  dat.s <- na.omit(dat.mf[,nf + 1:ns, drop=FALSE])

  if(is.null(upper)) upper <- ceiling(log(nrow(dat.s)) ^ 1.5)
  k.min <- ceiling(n*p / nf)

  stopifnot(upper >= k.min)

  x.t.hat <- Yminus(dat.mf, p-1)[1:nrow(dat.mf),]
  AIC <- rep(0,upper - k.min)

  for (k in k.min:upper){
    Y.minus <-  Yminus (dat.f, k)[1:nrow(dat.mf),]

    resids <- matrix(lm(dat.s ~ 0 + Y.minus[grid, ])$residuals, ncol=ns)
    len <- nrow(resids)

    num.par <- ns * nf * k
    Sigma_ML <- log(abs(det(crossprod(resids)/ len)))

    AIC[k - k.min + 1] <- Sigma_ML+  2 * num.par / len

  }
  k <- which(min(AIC) == AIC) + k.min - 1

  list(k=k, AIC=AIC)
}
