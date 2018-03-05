##################################################
##
##  MIXED FREQUENCY ESTIMATORS
##
## -----------------------------------------------
##
## 1. Estimation wrapper function
## 2. XYW Estimator
## 3. IVL Estimator
## 4. EM Algorithm for state space models
##
##################################################



#' Estimate AR parameters from MF-data
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

  ## AB: 2018-03-02
  if(method == 'XYW'){
    k <- args[['k']]
    a.hat.0 <- XYW(dat.mf, p=p, k=k)$a.hat
  }

  ## LK
  if(method == 'IVL'){
    k <- args[['k']]
    upper <- args[['upper']]
    a.hat.0 <- IVL(dat.mf, p=p, k=k, upper=upper)$a.hat
  }

  ## AB: 2018-03-02
  ## Note: EM algorithm estimates both, a.hat and sigma
  if (method == 'EM') {
      res.em <- EM (dat.mf, p, n.init=5, iter.max=100, tol=1e-4)
      a.hat.0 <- res.em$EM$a.hat
  }

  a.hat <- A.make.stable(a.hat.0, ...)$a

  sigma.help <- sigma.hat(dat.mf=dat.mf, a=a.hat, ...)


  list(a=a.hat, a.hat.0=a.hat.0, sigma=sigma.help$sigma, sigma.0=sigma.help$sigma.0)

}


#' XYW - Extended Yule Walker estimator for MF-data
#'
#' An estimator for the system-parameters of a VAR(p) process based on Extended-Yule Walker equations
#'
#' @param dat.mf The mixed-frequency data (T*n matrix), where the first nf components are the fast ones and the last ns components are the slow.
#' @param p Order of the AR polynomial
#' @param k A predefined k, which specifies the number columsn of \eqn{Z_0}. (TODO) If it is NULL, it will be estimated chosen such that the condition number of \eqn{Z_0} is minimal (the maximal number of columns of \eqn{Z_0}, denote by kmax, kmax = some number...)
#'
#' @return a.hat The estimator of the AR polynomial
#'
#' @export
XYW <- function(dat.mf, p, k=NULL ) {

    n <- ncol(dat.mf)
    detection <- detect_fast_slow(dat.mf)

    nf <- sum(detection$is.fast)
    ns <- n-nf
    N <- detection$N

    if(is.null(k)) {
        k <- n*p
    }
    
    autocovs <- stats::acf(dat.mf, na.action=na.pass, type="covariance", lag.max=k, demean=TRUE, plot=FALSE)$acf
      
    Z.1 <- XYW.Z1(autocovs, k, 1:nf)
    Z.0 <- XYW.Z0(autocovs, k, p, 1:nf)

    tryCatch({
        Aest <- Z.1%*%MASS::ginv(Z.0)
    }, error = function (e) print(e))
    
    list(a.hat=Aest)  
}

XYW.Z0 <- function (gam, k, p, fast) {

  Z <- NULL
  
  for (i in 0:(p-1)) {
    ith.row <- NULL
    
    for (j in 0:(k-1)) {
      if ((j-i)<0) if(length(fast)>1) gamma.j_i <- t(gam[abs(j-i)+1,fast,])
                   else gamma.j_i <- matrix(gam[abs(j-i)+1,fast,])
      else gamma.j_i <- gam[abs(j-i)+1,,fast]
      ith.row <- cbind(ith.row, gamma.j_i)
    }
    Z <- rbind(Z, ith.row)
  }
Z
}

XYW.Z1 <- function (gam,k,fast) {

  G.np <- NULL
  for(h in 1:k) { G.np <- cbind(G.np, gam[h+1,,fast]) }

  G.np
}




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

  if(is.null(k)) k <- AIC_IVL(dat.mf=dat.mf, p=p, n=n, nf=nf, N=N, grid=grid, upper=upper)$k

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
AIC_IVL <- function(dat.mf, p, n, nf, N, grid, upper=NULL){

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


#' The EM algorithm for state space models
#'
#' @param dat.mf The mixed-frequency data (T*n matrix), where the first nf components are the fast ones and the last ns components are the slow.
#' @param p Order of the AR polynomial
#' @param n.init Number of initial values, randomly generated
#' @param iter.max Maximum number of iterations
#' @param tol Convergence tolerance for algorithm
#'
#' @return MLE parameters a.hat and b.hat
#' 
#' @export
EM <- function (dat.mf, p, n.init, iter.max, tol=1e-4) {

    n <- ncol(dat.mf)
    detection <- detect_fast_slow(dat.mf)
    
    nf <- sum(detection$is.fast)
    ns <- n-nf
    N <- detection$N

    reml <- list()
    ll.names <-  NULL

    ## + length(thetas)
    loglik <- emtime <- rep(NA, 1+n.init)
    
    if (n.init>0) {
        for (j in 1:n.init) {
            ##
            ## hier: modgen(n, p)
            ##            rtau <- random_theta(m)$theta
            ##          tmp.ab <- theta.to.ab(rtau, m)
            tmp.ab <- modgen(n, p)
            
            ##print(tmp.ab)
            aborthere <- tryCatch({
                m.ss <- ar.to.ss(tmp.ab$a, tmp.ab$b, 0)
                FALSE
            }, error=function(e) { print(e); print(tmp.ab); TRUE })
                   
            ##print(aborthere)
            st <- NA
            
            if (!aborthere) { ## because of strange random starting value
                ##print(paste('EM for starting value', j))
                st <- system.time(reml[[j]] <- EM.algorithm(dat.mf, tmp.ab$a, tmp.ab$b, nf, N, iter.max, tol))[3]
            }
            if (!is.na(st)) {
                emtime[j] <- st
                loglik[j] <- reml[[j]]$loglik
            }
        }
        ll.names <- paste('r', 1:n.init, sep='')
    }
    
    ##if (!is.null(thetas) && length(thetas)>0) {
    ##    for (j in 1:length(thetas)) {
    ##        st <- NA
    ##        st <- system.time( reml[[k+j]] <- EM.algorithm(dat.mf, tmp.ab$a, tmp.ab$b, nf, N, iter.max, tol))[3]
    ##        if (!is.na(st)) {
    ##            emtime[k+j] <- st
    ##            loglik[k+j] <- reml[[k+j]]$loglik
    ##        }
    ##    }
    ##    ll.names <- c(ll.names, paste('s', 1:length(thetas), sep=''))
    ## }
    
    ##st <- NA
    ##st <- system.time (reml[[k+1+length(thetas)]] <- EM.algorithm(y, m, iter.max, tol, tau0=tau0))[3]
    ##  print(st)
    ##  print(reml[[k+1+length(thetas)]]$loglik)
    ##if (!is.na(st)) {
    ##    emtime[k+1+length(thetas)] <- st
    ##    loglik[k+1+length(thetas)] <- reml[[k+1+length(thetas)]]$loglik
    ## }
    ## ll.names <- c(ll.names, 'eyw')
    
                                        # print(loglik)
                                        # print(ll.names)
    names(loglik) <- ll.names
    names(emtime) <- ll.names
    j.max <- which.min(loglik)
                                        #print(j.max)
    ret <- reml[[j.max]]
    attr(ret, 'logliks') <- loglik
    
    list(EM=ret,
         reml=reml,
         loglik=loglik,
         emtime=emtime)
    
}



#' Fast Rcpp-EM Algorithm (EM0002)
#'
#' @export
#' @rdname EM
#'
EM.algorithm <- function (y, a, b, nf, N, iter.max, tol) {
  
    T <- nrow(y)
    n <- ncol(y)
    p <- floor(ncol(a)/n)
    q <- ncol(b)


    
    ## theta.to.ab - geht das auch anders?
    A <- A0 <- a
    b <- B0 <- b
    Sigma <- b %*% t(b)
  
    zeros <- matrix(rep(0,n*(p-1)*q), ncol=q)
    B <- rbind(b, zeros)
  
    ## STATE SPACE PARAMETER
    ## a bad random starting point could lead to P1-computation error:
    ## in this case a error must be thrown.
    m.ss <- ar.to.ss(A0, B0, 0)
    ##  csb <- parCtStB(m, dat$mf2)
  
    y[is.na(y)] <- 0
  
    ## LOGGING
    hist <- matrix(c(NA, NA, all(abs(eigen(A.companionform(A0))$values)<1), eigen(Sigma)$values[length(eigen(Sigma)$values)]), nrow=1)
    
    logLik1 <- logLik2 <- Inf
    ## logLikR1 <- logLikR2 <- Inf
  
    for (k in 1:iter.max) {
  
        res <- Estep0002(y, A.companionform(A), B, A, b%*%t(b), t(b), matrix(0,1, n*p), m.ss$V1, N, nf)

        hist[dim(hist)[1],1] <- res$loglik
        hist[dim(hist)[1],2] <- res$loglik<logLik2
    
        if (is.nan(res$loglik)) {
            cat("\n Error in loglik, NaN \n\n")
            print(k)
            print(b)
            return(list(A=A, B=B, b=b, a1=matrix(0,1,n*p), V1=m.ss$V1, loglik=res$loglik))
        }
        
        if (k>1 && emConverged(logLik2, res$loglik, tol)) {
            ## print ("converged")
            break
        }
        logLik2 <- res$loglik
        A <- res$Syx %*% solve(res$Sx)
    
        if (any(abs(eigen(A.companionform(A))$values)>=1)) {
            print("A unstable")
            print(abs(eigen(A.companionform(A))$values))
            
            print("10e-1")
            stab <- A.make.stable(a.unst=A,  makestable.tol=1e-1)
            A <- stab$a
            print("eigenvalues of stabilized A")
            print(abs(eigen(A.companionform(A))$values))
        }
        
        Sigma <- (1/(T-p))*(res$Sy - res$Syx %*% solve(res$Sx) %*% t(res$Syx))

        Sei <- eigen(Sigma)
        s <- q ##sum(Sei$values>1e-7)
        Se <- Sei$values[1:s]
        Se[Se<1e-8] <- 0
        b <- Sei$vector[,1:s] %*% sqrt(diag(Se, s))
    
        B <- rbind(b, zeros)
        
        m.ss <- ar.to.ss(a = A, b = b)
        
        
        tmp.out <- cbind(NA, NA, all(abs(eigen(A.companionform(A))$values)<1), eigen(Sigma)$values[length(eigen(Sigma)$values)])
        hist <- rbind(hist, tmp.out)
}
  
  colnames(hist) <- c('-logLik', 'll inc.', 'Astable', 'Sigma_1')
  rownames(hist) <- 1:dim(hist)[1]
  
  list(log=hist,
       a.hat=A,
       Sigma.hat=Sigma,
       iterations=k,
       loglik=res$loglik)
  
}


#' Has EM algorithm converged?
#'
emConverged <- function (ll0, ll1, threshold=1e-4) {

    avg <- (abs(ll0)+abs(ll1))/2

    d <- abs(ll1-ll0)

    (d/avg)<threshold

}


