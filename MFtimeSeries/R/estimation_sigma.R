## --------------------------------------------------
#' Estimation and projection of the innovation covariance matrix \eqn{\Sigma}
#'
#' Estimation of the innovation covariance matrix \eqn{\Sigma}  and
#' projection onto the space of positive-semidefinite matrices.
#' The estimation procedure and the projection is described in
#' L. Koelbl, A. Braumann, E. Felsenstein, M. Deistler (2016): Estimation of VAR Systems from Mixed Frequency Data: The Stock and the Flow Case, Advances in Econometrics, 35, pp. 43-73, DOI: 10.1108/S0731905320150000035002
#'
#' If q is not desired, it will be set automatically to n, i.e. dimension of the innovations. Only the stock case is supported up to now.
#' @param dat.mf The mixed-frequency data (T*n matrix)
#' @param a The estimated AR coefficients (n*np matrix)
#' @param q The desired rank of the innovation matrix \eqn{\Sigma} (n>=q>0)
#'
#' @return sigma The projected estimator of the innovation matrix \eqn{\Sigma}
#' @return sigma.0 The non-projected estimator of the innovation matrix \eqn{\Sigma}
#'
#' @export
sigma.hat <- function (dat.mf, a, ...) {
  n <- nrow(a)
  pars <- list(...)
  q <- pars[['q']]

  if(is.null(q)) q <- n

  stopifnot(q>0, n>=q)

  sigma.0 <- est.sigma(dat.mf, a)

  sigma <- proj.sigma(sigma.0, q)

  list(sigma = sigma, sigma.0 = sigma.0)
}

## --------------------------------------------------
#' Estimation of the innovation covariance matrix \eqn{\Sigma}
#'
#' Given the AR coefficients and the data,
#' an estimate for the innovation covariance matrix \eqn{\Sigma} is computed.
#'
#' @param dat.mf The mixed-frequency data (T*n matrix)
#' @param a The estimated AR coefficients (n*np matrix)
#'
#' @return sigma The non-projected estimator of the innovation matrix \eqn{\Sigma}
#' @seealso sigma.hat
#' @export
est.sigma <- function (dat.mf, a) {

  n <- nrow(a)
  p <- ncol(a) / n

  gamma0.hat <- calc.acf(dat.mf)[1,,]

  A <- A.companionform(a)
  G <- rbind(diag(n), matrix(0, ncol=n, nrow=n*(p-1)))

  s.vec <- solve(kronecker(t(G),t(G)) %*%
                   solve(diag((n*p)^2)-kronecker(A,A)) %*%
                   kronecker(G,G)) %*% matrix(gamma0.hat)

  sigma <- matrix(s.vec, ncol=n)

  sigma

}


