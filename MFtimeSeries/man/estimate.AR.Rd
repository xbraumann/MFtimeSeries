% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimator.R
\name{estimate.AR}
\alias{estimate.AR}
\title{Wrapper function for the estimation of stable AR and positive semi-definite innovation matrix}
\usage{
estimate.AR(dat.mf, p, method = "IVL", ...)
}
\arguments{
\item{dat.mf}{T*n data matrix of mf-data}

\item{p}{Order of the AR polynomial}

\item{method}{The method for the initialization}

\item{k}{A predefined k}

\item{upper}{The upper bound for the k-search}
}
\value{

}
\description{
Wrapper function for the estimation of stable AR and positive semi-definite innovation matrix
}
