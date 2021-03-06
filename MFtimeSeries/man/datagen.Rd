% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_data.R
\name{datagen}
\alias{datagen}
\title{Generate mixed-frequency data}
\usage{
datagen(T, a, b, N, nf, ...)
}
\arguments{
\item{T}{number of observations}

\item{a}{n*np matrix of system parameters}

\item{b}{n*q matrix of noise parameters}

\item{N}{mixed-frequency}

\item{nf}{number of slow components, ns := n - nf}

\item{...}{optional arguments, like burnin (standard = 2000)}
}
\description{
Aus Google-Dokument:
Generate mixed-frequency data, where DGP = AR(p) Prozess
}
\details{
Name: datagen
Parameter:
T, a, b (q = rank(b)), N, nf,
…
 burnin
 seed
 distribution of innovations (standard = rnorm)
 nf.index (which components are fast, …)?
 flow scheme
datagen

Generate mixed-frequency autoregressive time series data
}
