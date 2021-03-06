# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

glogdet <- function(X) {
    .Call('_MFtimeSeries_glogdet', PACKAGE = 'MFtimeSeries', X)
}

showMatrix <- function(X, name) {
    invisible(.Call('_MFtimeSeries_showMatrix', PACKAGE = 'MFtimeSeries', X, name))
}

#' The Kalman filter equations for state space systems
#' 
#' @param Y The data matrix
#' @param A The transiation matrix
#' @param B The matrix corresponding to the transition error, i.e. \eqn{\Sigma=BB'}
#' @param C The observation matrix
#' @param H The observation error covariance matrix
#' @param S The covariance between observation and transition error
#' @param a1 The mean of state x1
#' @param P1 The covariance of state x1
#' @param N The sampling frequency of the slow component
#' @param nf The number of fast components
#' 
#' @return A list with the usual Kalman filter variables.
#' 
kfilter <- function(Y, A, B0, C, H, S, a1, P1, N, nf) {
    .Call('_MFtimeSeries_kfilter', PACKAGE = 'MFtimeSeries', Y, A, B0, C, H, S, a1, P1, N, nf)
}

#' The Expectation step in the EM Algorithm
#' 
#' Compute the sufficient statistics in the EM Algorithm for state space models. Ignore the starting values.
#' 
#'
#' 
Estep0002 <- function(Y, A, B, C, H, S, a1, P1, N, nf) {
    .Call('_MFtimeSeries_Estep0002', PACKAGE = 'MFtimeSeries', Y, A, B, C, H, S, a1, P1, N, nf)
}

#' The Expectation step in the EM Algorithm
#' 
#' Compute the sufficient statistics in the EM Algorithm for state space models. Including starting values x1 with mean mu1 and covariance V1. See Shumway and Stoffer 1982.
#' 
#' 
#'
Estep0003 <- function(Y, A, B, C, H, S, a1, P1, N, nf) {
    .Call('_MFtimeSeries_Estep0003', PACKAGE = 'MFtimeSeries', Y, A, B, C, H, S, a1, P1, N, nf)
}

