% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation.R
\name{MakeMfdProcess}
\alias{MakeMfdProcess}
\title{Simulate a manifold-valued process}
\usage{
MakeMfdProcess(
  mfd,
  n,
  mu,
  pts = seq(0, 1, length.out = ncol(mu)),
  K = 2,
  lambda = rep(1, K),
  sigma2 = 0,
  xiFun = stats::rnorm,
  epsFun = stats::rnorm,
  epsBase = c("mu", "X"),
  basisType = c("cos", "sin", "fourier", "legendre01")
)
}
\arguments{
\item{mfd}{An object whose class specifies the manifold to use}

\item{n}{The number of curves}

\item{mu}{A matrix that specifies the mean function. Each column stands for a manifold-valued data point}

\item{pts}{The time points}

\item{K}{The number of components to simulate}

\item{lambda}{The eigenvalues}

\item{sigma2}{The variance for the isotropic measurement errors}

\item{xiFun}{The distribution of the RFPC scores}

\item{epsFun}{The distribution of the measurement errors}

\item{epsBase}{The base point of the (tangent) measurement errors.}

\item{basisType}{The type of basis on the tangent spaces}
}
\description{
Simulate a manifold-valued process
}
