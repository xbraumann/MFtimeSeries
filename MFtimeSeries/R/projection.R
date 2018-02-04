
#' Projection of a symmetric matrix on the space of positive semidefinte matrices.
#'
#' @param Sigma A matrix of dimension \eqn{n*n}{n\times n}.
#' @param q  q is the desired rank of sigma.
#' @param epsilon Negative eigenvalues will be set to epsilon.
proj.sigma<-function(sigma, q=NULL, epsilon=1e-7){

  n <- nrow(sigma)

  if(is.null(q) || q>n || q<1) q <- n

  eig <- eigen(sigma)
  eig.val <- cbind(eig$values[1:q], epsilon)
  eig.val <- apply(eig.val, 1, max)

  eig.vec <- matrix(eig$vectors[,1:q], nrow=n)

  sigma.proj <- eig.vec %*% diag(eig.val,nrow=q,ncol=q)  %*% t(eig.vec)

  sigma.proj

}

