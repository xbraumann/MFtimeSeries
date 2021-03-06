require(Matrix)


#' Generate randomly a model \eqn{a(z)}, b, sigma
#'
#' @param n Dimension of the process.
#' @param p Order of the autoregressive polynomial \eqn{a(z)}.
#'
#' @return a The AR polynomial in matrix form.
#' @return b
#' @return Sigma The variance matrix of the innovations.
#' @export
modgen <- function (n, p, ...) {

  stopifnot(n>1, p>0)

  a <- modgen.a(n, p, ...)

  b <- modgen.b(n, ...)

  sigma <- b%*%t(b)

  list(a=a, b=b, sigma = sigma)
}


## --------------------------------------------------
#' Randomly generate an AR polynomial \eqn{a(z)} in matrix form
#'
#' @param p The order of the autoregressive system.
#' @param n The dimension of the process.
#' @param complex.roots The number of complex roots of the polynomial \eqn{a(z)}.
#' @seealso gen_random_numbers
#'
#' @return a
#'
#'

modgen.a <- function (n, p, ...) {

  pars <- list(...)
  stopifnot(n>1, p>0)

  nb.c.roots <- pars[['complex.roots']]
  if (is.null(nb.c.roots)) nb.c.roots <- sample(seq(0, n*p, by=2), 1)
  else{
    stopifnot(n*p >= nb.c.roots, nb.c.roots >= 0,  nb.c.roots %% 2 == 0)
  }

  roots <- gen_random_numbers (n*p - nb.c.roots, nb.c.roots/2, ...)

  Lambda <- diag(x=roots$real,nrow=length(roots$real))

  if(nb.c.roots>0){
    complex.block <- lapply(roots$complex, function(x)matrix(c(Re(x),-Im(x),Im(x),Re(x)),2,2))

    Lambda <- as.matrix(Matrix::bdiag(Matrix::bdiag(complex.block),Lambda))
  }

  Trans <- matrix(0,ncol=n*p,nrow=n*p)
  Trans[1:n,] <- matrix(rnorm(n^2*p),nrow=n,ncol=n*p)

  if(p>1){
    Lambda.inv<- solve(Lambda)
    for(i in 1:(p-1)) {Trans[1:n+i*n,] <- Trans[1:n+(i-1)*n,] %*% Lambda.inv}
  }

  a <- Trans[1:n,] %*% Lambda %*% solve(Trans)

  a
}



#' Generate random numbers
#' @param r.roots The number of complex numbers, which should be generated
#' @param c.roots The number of complex numbers, which should be generated
#' @param real.r A vector of ranges of the moduli for the real numbers
#' @param complex.r A vector of ranges of the moduli for the complex numbers
#' @param ran A vector indicating the range, used if r.r an c.r are not specified
#'
#' @return Random real and complex random numbers
gen_random_numbers <- function (r.roots=2, c.roots=0, ...) {

  stopifnot(r.roots>=0, c.roots>=0)
  pars <- list(...)

  r.r <- pars[['real.r']]
  c.r <- pars[['complex.r']]

  if (length(c.r) != c.roots | length(r.r) != r.roots) {
    ran <- pars[['ran']]
    if (is.null(ran) | length(ran) != 2) ran <- c(0.7, 0.98)

    c.r <- runif(c.roots, ran[1], ran[2])
    r.r <- runif(r.roots, ran[1], ran[2])
  }

  c.phi <- runif(c.roots, -pi, pi)
  r.phi <- sample(c(-1,1), r.roots, replace=TRUE)

  c.z0 <- c.r * cos(c.phi) + sin(c.phi) * 1i
  r.z0 <- r.r * r.phi

  list(complex = c.z0, real = r.z0)
}



#' Generate a n*q matrix where columns have length 1
#'
#'
#' @param n The dimension of the process.
#' @param q The rank of sigma matrix.
#'
#' @return Matrix b of dimension n*q
modgen.b <- function (n, ...) {

  pars <- list(...)
  q <- pars[['q']]

  if (is.null(q)) q <- n
  stopifnot(n>=q, q>0)

  b <- matrix(rnorm(n*q), ncol=q, nrow=n)
  b <- apply(b,2,function(x){x/norm(x,type="2")})

  b
}
