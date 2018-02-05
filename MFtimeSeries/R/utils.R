

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
