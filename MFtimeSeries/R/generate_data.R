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
#'
#' @export
#' 
datagen <- function(T, a, b, N, nf, ...) {




}
