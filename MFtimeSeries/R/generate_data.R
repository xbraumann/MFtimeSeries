#'
#' Generate mixed-frequency data
#'
#' Aus Google-Dokument:
#' Generate mixed-frequency data, where DGP = AR(p) Prozess
#'
#' Name: datagen
#' Parameter:
#' T, a, b (q = rank(b)), N, nf,
#' …
#'  burnin
#'  seed
#'  distribution of innovations (standard = rnorm)
#'  nf.index (which components are fast, …)?
#'  flow scheme



#' datagen
#'
#' Generate mixed-frequency autoregressive time series data
#'
#' @param T number of observations
#' @param a n*np matrix of system parameters
#' @param b n*q matrix of noise parameters
#' @param N mixed-frequency
#' @param nf number of slow components, ns := n - nf
#' @param ... optional arguments, like burnin (standard = 2000)
#'
#' @export
#'
datagen <- function(T, a, b, N, nf, ...) {

    require(dse, quietly = TRUE) ## we use ARMA, simulate
    n <- nrow(a)
    stopifnot(N>=1, n==nrow(b), nf>=1, n>nf)
    ## stopifnot(qr(b)$rank==ncol(b))

    pars <- list(...)

    if (!is.null(pars[['burnin']])) T.start <- pars[['burnin']]
    else T.start <- 2000

                                        # burn-in
    T1 <- T+T.start

    n <- nrow(a)
    q <- ncol(b)  ## or qr(b)$rank

    ## create ARMA object for the simulate function
    AR <- array.reorder((-1)*matrix.to.array(a), a0=FALSE)
    ## NEW VERSION with base::aparm: test!
    ## AR <- aperm(matrix.to.array(cbind(diag(1, n), -a)), c(3, 1, 2))
    MA <- diag(1, n)
    arma.system <- ARMA(A=AR, B=MA)

    ## innovations - normal i.i.d
    U <- matrix(rnorm(T1*q, 0, 1), ncol=q, byrow=TRUE)
    V <- U %*% t(b)

    data <- simulate( arma.system, sampleT=T1, noise=V )
    ## the T*n data-matrix
    Y.0 <- data$output[-(1:T.start),,drop=FALSE]

    ## mixed frequency data
    y  <- Y.0
    y[c(rep(TRUE,N-1),F) , (nf+1):n] <- NA

    ## RETURN:
    ret <- list(high=Y.0,
                mixed=y)

    ret

}

## array.reorder DELETE?!
#' is base::aperm doing the same?
#' array.reorder(matrix.to.array(cbind(diag(1, 3), -a)), a0=TRUE)-aperm(matrix.to.array(cbind(diag(1, 3), -a)), c(3, 1, 2))
#'
#' convert a 3-dim array with dims [n,n,lag] into array with dims [lag,n,n]
#' a0=FALSE: add a0=I_n
#' only used to convert m$a into input of dse::ARMA
#'
#' @export
array.reorder <- function (a, a0=FALSE, append.a0.only=FALSE) {

  s <- 1*a0

  if (!s) a1 <- array(c(diag(1, dim(a)[1]), a), dim=c(dim(a)[1], dim(a)[2], dim(a)[3]+1))
  else a1 <- a

  if (append.a0.only) return (a)

  b <- array(rep(0, dim(a)[1]*dim(a)[2]*(dim(a)[3]+(1-s))), dim=c(dim(a)[3]+(1-s), dim(a)[1], dim(a)[2]))

  for (x in 1:dim(a1)[1]) {

    for (y in 1:dim(a1)[2]) {

      for (z in 1:dim(a1)[3]) b[z, x, y] <- a1[x, y, z]

    }
  }

  b
}


#' convert matrix to array
#' @export
matrix.to.array <- function (a) {

   p <- dim(a)[2]/dim(a)[1]
   n <- dim(a)[1]

   A <- array(0, dim=c(n,n,p))

   for (i in 1:p) {A[,,i] <- a[,(1+(i-1)*n):(i*n)] }

   A
}
