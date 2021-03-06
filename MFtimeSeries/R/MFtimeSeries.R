#' VAR estimation for mixed frequency data
#'
#' This package provides methods which allow to estimate the system and noise parameters of a
#' multivariate auto-regressive system in case of mixed frequency data.
#'
#' We assume that the data generating system is given by
#' \deqn{y_t = A_1 y_{t-1} + \ldots + A_p y_{t-p} + b \varepsilon_t}
#' where \eqn{\varepsilon_t} is white noise with mean 0 and covariance \eqn{I_q}, \eqn{b} is a
#' matrix of dimension \eqn{n*q}{n\times q}, and \eqn{A_i} are matrices of dimension \eqn{n*n}{n\times n}.
#'
#' The package provides the following:
#' \describe{
#'  \item{Estimation methods}{The XYW estimator corresponding to the pseudoinverse of \eqn{Z_0}
#'  (see \code{\link{EYW}}), the GMM estimator (see \code{\link{GMM}}), the Gaussian MLE based on
#'  the EM algorithm (see \code{\link{EM}}).}
#'  \item{Stabilization of the \eqn{A} parameters}{\code{\link{A.make.stable}}, \code{\link{A.make.stable.convex.approx}}.}
#'  \item{Asymptotic covariances of the XYW estimator}{see \code{\link{Bartlett}}.}
#'  \item{Mixed frequency forecasting}{based on the Kalman filter.}
#' }
#'
#' @name MFtimeSeries-package
#' @aliases MFtimeSeries-package MFtimeSeries
#' @docType package
#' @useDynLib MFtimeSeries
#' @importFrom Rcpp sourceCpp
NULL
