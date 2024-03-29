% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NoiseTangent.R
\name{NoiseTangent}
\alias{NoiseTangent}
\alias{NoiseTangent.Sphere}
\alias{NoiseTangent.Euclidean}
\alias{NoiseTangent.SPD}
\title{Generate noise on the tangent space}
\usage{
NoiseTangent(mfd, n, mu, sigma2 = 0, epsFun = rnorm)

\method{NoiseTangent}{Sphere}(mfd, n, mu, sigma2 = 0, epsFun = stats::rnorm)

\method{NoiseTangent}{Euclidean}(mfd, n, mu, sigma2 = 0, epsFun = stats::rnorm)

\method{NoiseTangent}{SPD}(mfd, n, mu, sigma2 = 0, epsFun = stats::rnorm)
}
\arguments{
\item{mfd}{A manifold object created by \code{\link{createM}}}

\item{n}{Sample size}

\item{mu}{A matrix or a vector containing the base point(s)}

\item{sigma2}{Noise variance}

\item{epsFun}{The distribution for generating the univariate errors with unit variance for the tangent space coordinates}
}
\value{
An array with size (n, dimIntrinsic, nPts)
}
\description{
For each time point (column in mu), generate n isotropic noises distributed as epsFun on the tangent space at mu(t)
}
\section{Methods (by class)}{
\itemize{
\item \code{NoiseTangent(Sphere)}: Method

\item \code{NoiseTangent(Euclidean)}: Method

\item \code{NoiseTangent(SPD)}: Method

}}
