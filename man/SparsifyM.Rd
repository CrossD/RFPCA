% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation.R
\name{SparsifyM}
\alias{SparsifyM}
\title{Sparsify a dense Riemannian process}
\usage{
SparsifyM(X, pts, sparsity)
}
\arguments{
\item{X}{A three dimensional array such that X[i, j, t] = X_{ij}(t).}

\item{pts}{A vector of grid points corresponding to the columns of \code{samp}.}

\item{sparsity}{A vector of integers. The number of observation per sample is chosen to be one of the elements in sparsity with equal chance}
}
\description{
Sparsify a dense Riemannian process
}
