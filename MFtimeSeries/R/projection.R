#' AR(p) nearest stable system using succesive convex approximation
#'
#' This algorithm has been described in Koelbl et.al. (2016): Estimation of VAR Systems from Mixed Frequency Data: The Stock and the Flow Case (Advances in Econometrics) and uses the idea of
#' Orbandexivry et.al. (2013). Nearest stable system using successive convex approximations (Automatica).
#' The following initialization methods are available:
#' \describe{
#'   \item{EV.ref1}{Uses the idea of Blaschke for reflecting the unstable roots of the AR polynomial}
#'   \item{EV.ref2}{Only reflects the unstable roots of the AR polynomial}
#'   \item{EV.scale}{Scales the unstable roots of the AR polynomial to \eqn{1-epsilon}
#'   \item{A.scale}{Scales the AR parameters directly}
#' }
#'
#'
#' @param a.unstable The unstable AR parameters in matrix form (n*np)
#' @param init.method The different initialization methods 'EV.ref1', 'EV.ref2', 'EV.scale', 'A.scale'
#' @param alpha The size of the Dikin Elipsoid
#' @param epsilon Convergence criterion
#' @param makestable.tol The minimal length of step $h$
#'
#' @return A list with the following elements
#' \describe{
#'   \item{a}{The nearest stable AR polynomial}
#'   \item{init}{The AR polynomial used for the initialization}
#'   \item{lam}{The step size at each iteration}
#'   \item{norm.step}{The norm of the matrix which has been added at each iteration step}
#'   \item{steps}{The number of steps needed until a stable \eqn{A} has been reached}
#' }
#'
#' @export
A.make.stable <- function (a.unstable, ...){

  require(matrixcalc, quietly = TRUE)

  if(A.stable(a.unstable)) return(list(a = a.unstable, init = a.unstable, lam=0,
                                   norm.step=0, steps=0))

  pars <- list(...)
  n <- nrow(a.unstable)
  p <- ncol(a.unstable) / n

  init.method <- pars[['init.method']]
  if (is.null(init.method)) init.method <- "EV.ref1"

  iterations <- pars[['iterations']]
  if (is.null(iterations)) iterations <- 50

  tolerance <- pars[['tolerance']]
  if (is.null(tolerance)) tolerance <- 1e-4

  alpha <- pars[['alpha']]
  if (is.null(alpha)) alpha <- 0.9999

  init <- pars[['init']]
  if (is.null(init)) init <- initialization(a.unstable = a.unstable, init.method = init.method, ...)



  A.0 <- A.companionform(init)

  j <- 0
  step <- matrix(1)
  lam <- norm.step <- rep(0, iterations)

  #' as S is open, we have to specify when to stop by tolerance
  while(j < iterations && norm(step, type="F") > tolerance){

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

  list(a = A.0[1:n,], init = init, lam = lam,
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


#' Different initialization algorithms for the stabilization algorithm for the AR polynomial
#'
#' The following methods are available:
#'
#' \describe{
#'   \item{EV.ref1}{Uses the idea of Blaschke for reflecting the unstable roots of the AR polynomial}
#'   \item{EV.ref2}{Only reflects the unstable roots of the AR polynomial}
#'   \item{EV.scale}{Scales the unstable roots of the AR polynomial to \eqn{1-epsilon}
#'   \item{A.scale}{Scales the AR parameters directly}
#' }
#'
#' @param a.unstable The unstable AR system parameters in matrix form (n*np)
#' @param init.method The name of the initialization method 'EV.ref1', 'EV.ref2', 'EV.scale', 'A.scale'
#'
#' @return a.stable The calculated stable AR polynomial
#'
#' @export
#'
initialization <- function (a.unstable, init.method = 'EV.ref1', ...){

  if(A.stable(a.unstable)) return(a.unstable)

  pars <- list(...)
  n <- nrow(a.unstable)
  p <- ncol(a.unstable)/n

  if(init.method == 'EV.ref1') a.stable <- Blaschke.transformation(a.unstable)$a
  else{

    eig     <- eigen(A.companionform(a.unstable))
    eig.val <- eig$values
    eig.vec <- eig$vector

    epsilon <- pars[["init.epsilon"]]
    if(is.null(epsilon)) epsilon <- 1e-6

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

      for(i in 1:p){ a.stable[, 1:n + (i-1)*n] <- a.stable[, 1:n + (i-1)*n] * mu^i  }
    }
  }

  a.stable
}


#' Stabilize estimator of A using the idea of Blaschke.
#'
#' We are assuming that a(z) has different zeros. The idea is to reflect the roots which are outside the unit circle.
#'
#'
#' @param a.unstable The unstable AR system parameters in matrix form (n*np)
#' @param sigma The variance matrix \eqn{\Sigma} of the innovations.
#'
#' @return a The stable polynomial
#' @return sigma The corresponding sigma
#'
#' @export
#'
Blaschke.transformation <- function(a.unstable, sigma=NULL){

  n <- nrow(a.unstable)
  p <- ncol(a.unstable)/n
  A <- A.companionform(a.unstable)

  if(is.null(sigma)) sigma <- diag(n)

  b <- t(chol(sigma))

  a.bla <- solve(b) %*% cbind(diag(n), -a.unstable)


  while(!A.stable(A)){

    eig.val <- eigen(A)$values

    z.0 <- 1/eig.val[1]
    a.z0 <- a.bla[,1:n]

    for(i in 1:p){ a.z0 <- a.z0 + z.0^i * a.bla[, 1:n + n*i] }
    eig.vec <- solve( eigen(a.z0)$vectors )

    K <- rbind(eig.vec[n,], cbind(0, diag(n-1)))
    K <- t( orthonormalization_complex(t(K)) )

    a.bla <- K %*% a.bla

    # enumerator
    enumerator <- c(a.bla[1,], rep(0,n)) -c(rep(0,n), a.bla[1,]) * Conj(z.0)
    a.bla[1,] <- 0

    #denumerator
    for(j in 0:p){
      for(k in (j+1):(p+1)){
        a.bla[1, 1:n + n*j] <- a.bla[1, 1:n + n*j] + z.0^(k-j-1) * enumerator[1:n + n*k]
      }
    }


    A.0.inv <- solve(a.bla[, 1:n])

    A <- A.companionform(-A.0.inv %*% a.bla[, (n+1):ncol(a.bla)])
  }

  if( norm(Im(A[1:n,]), type="F") < 10^-10 ) A <- Re( A[1:n,] )
  else print("A is complex")

  b <- A.0.inv
  sigma <- b %*% t(Conj(b))

  if( norm(Im(sigma), type="F") < 10^-10 ) sigma <- Re(sigma)
  else print("Sigma is complex")


  list(a = A,sigma = sigma)

}

#' @describeIn A.make.stable Orthonomalization of a matrix
orthonormalization_complex <- function (u = NULL, basis = TRUE, norm = TRUE){
  if (is.null(u))
    return(NULL)
  if (!(is.matrix(u)))
    u <- as.matrix(u)
  p <- nrow(u)
  n <- ncol(u)
  if (prod(abs(La.svd(u)$d) > 1e-08) == 0)
    stop("colinears vectors")
  if (p < n) {
    warning("too much vectors to orthogonalize.")
    u <- as.matrix(u[, 1:p])
    n <- p
  }

  v <- u
  if (n > 1) {
    for (i in 2:n) {
      coef.proj <- c(crossprod(u[, i], v[, 1:(i - 1)]))/diag(crossprod(v[,1:(i - 1)],Conj(v[,1:(i - 1)])))
      v[, i] <- u[, i] - matrix(v[, 1:(i - 1)], nrow = p) %*%
        matrix(Conj(coef.proj), nrow = i - 1)
    }
  }

  if (norm) {
    coef.proj <- 1/sqrt(diag(crossprod(v,Conj(v))))
    v <- t(t(v) * coef.proj)
  }
  return(v)
}


#' Projection of a symmetric matrix on the space of positive semi-definite matrices.
#'
#'
#' @param sigma A symmetric matrix of dimension n*n.
#' @param q  q is the desired rank of sigma.
#' @param epsilon Negative eigenvalues will be set to epsilon.
#'
#' @return sigma.proj The projected estimator of the innovation matrix.
proj.sigma <- function(sigma, q=NULL, epsilon=1e-7){

  n <- nrow(sigma)

  if(is.null(q) || q>n || q<1) q <- n

  eig <- eigen(sigma)
  eig.val <- eig$values[1:q]
  eig.val[eig.val < epsilon]  <- epsilon

  eig.vec <- matrix(eig$vectors[,1:q], nrow=n)

  sigma.proj <- eig.vec %*% diag(eig.val,nrow=q,ncol=q)  %*% t(eig.vec)

  sigma.proj

}

