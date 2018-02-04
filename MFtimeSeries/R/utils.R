#' @rdname systemparameters
#' @export
A.companionform <- function(a) {
  stopifnot(is.matrix(a))

  n <- nrow(a)
  np <- ncol(a)

  A <- rbind(a, cbind(diag(np-n),matrix(0, nrow=(np-n), ncol=n)))
  A
}
