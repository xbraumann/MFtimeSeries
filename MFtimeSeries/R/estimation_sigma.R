## --------------------------------------------------
#' Estimation of the innovation covariance matrix and
#' projection onto the space of positive-semidefinite matrices
#'
#' @param data The mixed-frequency data
#' @param a The estimated AR coefficients
#' @param q The desired rank of the estimated innovation matrix
#'
#' @return sigma The projected estimate of sigma
#' @return sigma.0 The non-projected estimate of sigma
#'
#' @export
sigma.hat <- function (dat.mf, a, q=NULL, ...) {
  n <- nrow(a)
  if(is.null(q)) q <- n

  stopifnot(q>0, n>=q)

  sigma.0 <- est.sigma(dat.mf,a)

  sigma <- proj.sigma(sigma.0, q)

  list(sigma = sigma, sigma.0 = sigma.0)
}

## --------------------------------------------------
#' Compute the innovation covariance matrix
#'
#' Given the AR coefficients and the data,
#' an estimate for the innovation covariance matrix \eqn{\Sigma} is computed.
#'
#' @param data The mixed-frequency data
#' @param a The estimated AR coefficients
#'
#' @export
est.sigma <- function (data, a) {

  n <- nrow(a)
  p <- ncol(a) / n

  gamma0.hat <- calc.acf(data)[1,,]

  A <- A.companionform(a)
  G <- rbind(diag(n), matrix(0, ncol=n, nrow=n*(p-1)))

  s.vec <- solve(kronecker(t(G),t(G)) %*%
                   solve(diag((n*p)^2)-kronecker(A,A)) %*%
                   kronecker(G,G)) %*% matrix(gamma0.hat)

  sigma <- matrix(s.vec, ncol=n)

  sigma

}


