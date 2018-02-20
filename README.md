# MFtimeSeries: Mixed-Frequency Time Series

This is an R-package implementation based on the theory of a couple of research papers which deals with the problem of analyzing autoregressive systems in the case of mixed-frequency observations. For more details about the theory, see the References section below.

Have fun,
Alexander Braumann and Lukas Kölbl

## Introduction
This package deals with the estimation of multivariate autoregressive systems in the case of observations with different sampling rates, the so-called mixed-frequency (MF) observations. In this package the case where the output variable can be separated into a fast (high-frequency) and a slow (low-frequency) component will be given attention. 

An example of such an observation scheme would be a two-dimensional time series in which the first component is observed monthly, for instance unemployment, and the second quarterly, for instance GDP. It is assumed that the underlying system generates the output at each time point, the so-called high-frequency, however, the output of the slow component is only observed at an integer multiple of the high-frequency. Here two cases for the slow component are considered: the stock and the flow case.

This special observation pattern plays an important role in high-dimensional time series, where the availability of the univariate time series is given at different sampling rates. A popular approach to model high-dimensional time series are the generalized linear dynamic factor models, where the static factors can be typically modeled by a singular autoregressive system. 

## Usage:
For a quick demo see the demo.R file.

The following estimation procedures are supported:
* XYW: Extended Yule Walker Estimator
* IVL: Instrumental Variable Estimator
* EM: EM-Algorithm

Furthermore, a stabilization algorithm for the AR parameters as well as for the variance covariance matrix is supported.

## References:
* L. Koelbl, M. Deistler (2018): A New Approach for Estimating VAR Systems in the Mixed-Frequency Case, Statistical Papers, 10.1007/s00362-018-0985-1
* M. Deistler, L. Koelbl, B.D.O. Anderson (2017): Non-Identifiability of VMA and VARMA Systems in the Mixed Frequency Case, Econometrics and Statistics, 4, pp. 31-38, DOI: 10.1016/j.ecosta.2016.11.006
* B.D.O. Anderson, M. Deistler, E. Felsenstein, L. Koelbl (2016): The Structure of Multivariate AR and ARMA Systems: Regular and Singular Systems; the Single and the Mixed Frequency Case, Journal of Econometrics, 192, Issue 2, pp. 366-373, DOI: 10.1016/j.jeconom.2016.02.004
* L. Koelbl, A. Braumann, E. Felsenstein, M. Deistler (2016): Estimation of VAR Systems from Mixed Frequency Data: The Stock and the Flow Case, Advances in Econometrics, 35, pp. 43-73, DOI: 10.1108/S0731905320150000035002
* L. Koelbl (2015): VAR Systems: g-Identifiability and Asymptotic Properties of Parameter Estimates for the Mixed-Frequency Case, Dissertation, Technische Universität Wien
* B.D.O. Anderson, M. Deistler, E. Felsenstein, B. Funovits, L. Koelbl, M. Zamani (2015): Multivariate AR Systems and Mixed Frequency Data: g-Identifiability and Estimation, Econometric Theory, 32, Issue 4, pp. 793-826, DOI: 10.1017/S0266466615000043

## Planned updates/extensions:

This is a work in progress.  Planned extensions include:
* **Flow Case**
   * Extension of the algorithms to the mixed-frequency flow case.
   
* **ARMA**
   * Extension of the algorithms to the ARMA case. 

## Major updates and releases:
* 03/20/2018: Initial release, version 1.0.

