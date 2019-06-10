# Multivariate eigen analysis based on a smoothed covariance surface.
# smoothCov: A 4-dimensional array, whose dimensions corresponds to (t, s, j, l), where j, l are entry indices in the ambient space.
# K: The number of components to select. If not specified, default to all components
# FVE The desired FVE, out of the total specified by the first K components
#
# Returns:
# lambda: a vector of eigenvalues
# phi: a matrix of eigenvectors, where the jth column contains [phi_j1^T, ..., phi_jp^T] (the time indices are consecutive)

EigenAnalysis <- function(smoothCov, regGrid, K, FVE=1, muWork=NULL, tol=1e-14, verbose=TRUE, fastEig=FALSE) {

  if (missing(K)) {
    K <- Inf
  } else {
    stopifnot(K >= 1)
  }

  if (length(regGrid) == 1) { # For scalar HS/L2 cases
    gridSize <- 1
  } else {
    gridSize <- regGrid[2] - regGrid[1]
  }
  nT <- dim(smoothCov)[1]
  p <- dim(smoothCov)[3]
  
  mat <- matrix(aperm(smoothCov, c(1, 3, 2, 4)), nT * p, nT * p)
  mat <- (mat + t(mat)) / 2 # In case smoothCov is not symmetrical due to numerics
  if (fastEig) {
    eig <- tryCatch({
      RSpectra::eigs_sym(mat, min(K, nrow(mat)))
    }, error=function(e) {
      eigen(mat)
    })
      
  } else {
    eig <- eigen(mat)
  }

  positiveInd <- eig[['values']] / (tol + eig[['values']][1]) > tol

  if (sum(positiveInd) == 0) {
    stop('All eigenvalues are negative. The covariance estimate is incorrect.')
  }

  d <- eig[['values']][positiveInd]
  eigenV <- eig[['vectors']][, positiveInd, drop=FALSE]

  if (K > length(d)) {
    if (verbose && !is.infinite(K)) {
      warning(sprintf("Returning only %d < K = %d components with positive eigenvalues.\n", length(d), K)) 
    }
    K <- length(d)
  } 
    
  d <- d[seq_len(K)]
  eigenV <- eigenV[, seq_len(K), drop=FALSE]
  cumFVE <- cumsum(d) / sum(d)
  FVEK <- min(which(cumFVE >= FVE))
  d <- d[seq_len(FVEK)]
  eigenV <- eigenV[, seq_len(FVEK), drop=FALSE]

  # normalization
  if (is.null(muWork)) {
    muWork = seq_len(nrow(eigenV))
  }
  
  phi <- apply(eigenV, 2, function(x) {
                 x <- x / sqrt(gridSize)
                 if ( 0 <= crossprod(x, muWork) ) {
                   x
                 } else {
                   -x
                 }
  })
  lambda <- gridSize * d

  fittedCov <- phi %*% diag(x=lambda, nrow = length(lambda)) %*% t(phi)
  fittedCov <- aperm(array(fittedCov, c(nT, p, nT, p)), c(1, 3, 2, 4))

  return(list(lambda = lambda, phi = phi, 
              cumFVE = cumFVE, fittedCov=fittedCov))
}

