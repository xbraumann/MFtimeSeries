
#' state space representation of varp system
#'
#' @param R_ii The diagonal element of the measurement covariance matrix \eqn{R}
#' @param scheme Weights for flow data
#' 
#' @rdname systemparameters
#' @export
ar.to.ss <- function (a, b, R_ii=NULL, scheme=NULL, nf=NULL) {
 
    n <- dim(a)[1]
    p <- floor(dim(a)[2]/n)
    q <- ncol(b)

    if (is.null(scheme)) m <- p
    else m <- max(p, length(scheme))
    
    if (is.null(R_ii)) R <- diag(1e-10,n)
    else R <- diag(R_ii, n)
    
    Sigma <- b%*%t(b)

    S <- matrix(0, n*m, n*m)
    S[1:n, 1:n] <- Sigma

    if (is.null(scheme)) {
        oH <- matrix(rep(diag(n)*0, p-1), nrow=n)
        H <- cbind(diag(n), oH)
    } else {

        H <- cbind(diag(1, nf), matrix(0, nf, m*n-nf))
        H <- rbind(H, kronecker(matrix(scheme, nrow=1), cbind(matrix(0,n-nf, nf), diag(1,n-nf))))
    }

    A <- A.companionform(cbind(a, matrix(0, n, (m-p)*n)))
    
    P1 <- matrix(solve(diag(n^2*m^2)-kronecker(A, A))%*%matrix(S), ncol=n*m)
    mu1 <- as.vector(rep(0, n*m))

    zeros <- matrix(rep(0,n*(p-1)*q), ncol=q)
    
    list(A=A,
         B=rbind(b, zeros),
         C=H,
         a1=mu1,
         V1=P1,
         Q=S,
         R=R)
   
}




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
  sel <- min(nrow(x),500)
  grid <- which(apply(x[1:sel,], 1, function(x) {all(!is.na(x)) }))
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

