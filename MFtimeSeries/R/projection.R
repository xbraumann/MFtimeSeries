
#' Make AR paramter estimates stable
#'
#' ARp - Nearest stable system using succesive convex approximation
#'
#' The following initialization methods are available:
#' \describe{
#'  \item{ref}{Reflecting the roots, lieing inside the unit circle, on the unit circle.}
#'  \item{scale1}{...}
#'  \item{scale2}{...}
#'  \item{unit}{...}
#' }
#'
#'
#' @param a.unst The unstable AR coefficients
#' @param init.method Methods a)-d) (see Details)
#' @param alpha Has to do with the size of the Dikin Elipsoid
#' @param epsilon Convergence criterion
#' @param makestable.tol The minimal length of step $h$
#'
#' @return A list with the following elements
#' \describe{
#'   \item{a}{The nearest stable AR coefficients.}
#'   \item{init}{...}
#'   \item{lam}{...}
#'   \item{nor}{...}
#'   \item{nor.h}{...}
#'   \item{steps}{The number of steps needed until a stable \eqn{A} has been reached.}
#' }
#'
#' @export
A.make.stable <- function (a.unstable, init.method='EV.scale', alpha=0.9999, tolerance=1e-4, iterations=50, ...){

  if(A.stable(a.unstable)) return(list(a = a.unstable, init = a.unstable, lam=0,
                                   norm.step=0, steps=0))

  pars <- list(...)
  n <- nrow(a.unstable)
  p <- ncol(a.unstable) / n

  #' start with a stable polynomial
  init <- A.0 <- pars[['A.0']]

  if (is.null(A.0) | nrow(A.0) != ncol(A.0)) {
    init <- A.0 <- initialization(a.unstable = a.unstable, init.method = init.method, ...)
  }

  step <- j <-0
  lam <- norm.step <- rep(0, iterations)

  #' as S is open, we have to specify when to stop by tol
  while(j < iterations && norm(step, type="F") > tol){

    j <- j + 1
    #' compute the optimal P* and Q, (6.1.5-6)
    PQ.star <- optimal.PQ(A.0, n, p)
    Q.inv   <- PQ.star$Q.inv
    P       <- PQ.star$P

    #' Compute the operator matrix B (6.1.8)
    K <- matrixcalc::commutation.matrix(r=n, c=n*p)
    B <- 2*(kronecker(P + P %*% t(A.0) %*% Q.inv %*% A.0 %*% P, Q.inv[1:n,1:n]) +
            kronecker(P %*% t(A.0) %*% Q.inv[,1:n], Q.inv[1:n,] %*% A.0 %*% P) %*% K)

    B   <- (B + t(B)) / 2
    eig <- eigen(B)
    e   <- t(eig$vectors)
    d   <- eig$values

    enumerator <- d * (e %*% matrix(a.unstable - A.0[1:n,]) )^2

    #' optimize psi(lambda) (6.1.10)
    optout <- stats::uniroot(f = psi.fn, interval = c(0,1e100),
                            extendInt = "no", tol = 1e-10,
                            d = d, enumerator = enumerator, alpha = alpha)
    lambda <- optout$root

    step <- matrix( solve( diag(n^2*p) + lambda * B) %*% matrix(a.unstable - A.0[1:n,]), nrow=n)

    A.0[1:n,] <- A.0[1:n,] + step
    norm.step[j] <- norm(step, type="F")
    lam[j] <- lambda
  }

  norm.step <- norm.step[1:j]
  lam       <- lam[1:j]

  list(a = A.0[1:n,], init = init[1:n,], lam = lam,
       norm.step = norm.step, steps = j)

}

#' @describeIn A.make.stable function that calculates optimal P and Q
optimal.PQ <- function (A, n, p){

  Q.inv <- matrix(
    solve( diag((n*p)^2) - kronecker(t(A),t(A)) ) %*%
      matrix( n*p*diag(n*p) ), nrow = n*p)
  Q.inv <- (Q.inv + t(Q.inv)) / 2

  Q <- solve(Q.inv)
  P <- matrix(
    solve( diag((n*p)^2) - kronecker(A,A) ) %*%
      matrix(Q), nrow = n*p)

  P <- (P + t(P)) / 2
  list(P=P, Q.inv=Q.inv)
}

#' @describeIn A.make.stable objective function for stats::uniroot
psi.fn <- function (lambda, d, enumerator, alpha=0.999){
  sum(enumerator/ (1 + lambda * d)^2) - alpha
}


#' @describeIn A.make.stable Initialization of the algorithm
initialization <- function (a.unstable, init.method='EV.scale', ...){

  if(A.stable(a.unstable)) return(a.unstable)

  pars <- list(...)
  n <- nrow(a.unstable)
  p <- ncol(a.unstable)/n

  eig     <- eigen(A.companionform(a.unstable))
  eig.val <- eig$vector
  eig.vec <- eig$values

  epsilon <- pars[["init.epsilon"]]
  if(is.null(epsilon)) epsilon <- 1e-6

   if(init.method=='EV.ref1') {
     #OPEN
  #   a.stable<-A.make.stable(a.unst)$a
     a.stable<-a.unstable
  }

  if(init.method == 'EV.scale' || init.method == 'EV.ref2') {

    not.stable <- abs(eig.val)>1

    if(init.method == 'EV.scale') eig.val[not.stable]  <- (1 - epsilon) * exp(1i * Arg( eig.val[not.stable] ))
    else eig.val[not.stable]  <- 1 / Conj(eig.val[not.stable])

    Trans <- NULL
    for(i in 1:p){
      Trans <- rbind(Trans, eig$vector[1:n,] %*% diag( eig.val^(-i+1) ))
    }

    a.stable <- Re(Trans %*% diag(eig.val) %*% solve(Trans))[1:n,]
  }

  if(init.method == 'A.scale'){
    mu <- min(abs( 1 / eig.val)) * (1-epsilon)
    a.stable <- a.unstable

    for(i in 1:p){
      a.stable[, 1:n + (i-1)*n] <- a.stable[, 1:n + (i-1)*n] * mu^i
    }
  }

  a.stable
}


#' Projection of a symmetric matrix on the space of positive semidefinte matrices.
#'
#' @param Sigma A matrix of dimension \eqn{n*n}{n\times n}.
#' @param q  q is the desired rank of sigma.
#' @param epsilon Negative eigenvalues will be set to epsilon.
proj.sigma<-function(sigma, q=NULL, epsilon=1e-7){

  n <- nrow(sigma)

  if(is.null(q) || q>n || q<1) q <- n

  eig <- eigen(sigma)
  eig.val <- eig$values[1:q]
  eig.val[abs(eig.val) < epsilon]  <- epsilon

  eig.vec <- matrix(eig$vectors[,1:q], nrow=n)

  sigma.proj <- eig.vec %*% diag(eig.val,nrow=q,ncol=q)  %*% t(eig.vec)

  sigma.proj

}

