# Estimate multivariate FPC scores through conditional expectations.
#
# y: A list of function values
# t: A list of time points corresponding to y
# mu: A matrix of size p by nT
# obsGrid: supports fittedCov and phi. May be a truncated version.
# fittedCov: An array of size nT by nT by p by p
# lambda: A vector of eigenvalues
# phi: An array of size nT by p by K
# sigma2: A scalar or a vector of nonnegative scalars, corresponding to the process error variances (or regularization)
# Assumes for mu, phi, and fittedCov, each tVec in t are all supported on obsGrid.
# return: ret is a 3 by n array, with the first row containing the xiEst, second row containing the xiVar, and the third containing the fitted values. 
CEScores <- function(y, t, optns, mu, obsGrid, fittedCov, lambda, phi, sigma2) {

  if (length(lambda) != dim(phi)[3])
    stop('Number of eigenvalues is not the same as the number of eigenfunctions.')

  if (is.null(sigma2))
    sigma2 <- 0

  Sigma_Y <- fittedCov + 
    outer(diag(1, dim(fittedCov)[1]), 
          diag(sigma2, dim(fittedCov)[3]))

  MuPhiSig <- MuPhiSigMult(t, obsGrid, mu, phi, Sigma_Y)
  ret <- mapply(function(yVec, muphisig) {
                  # browser()
           fdapace:::GetIndCEScores(
             yVec, muphisig$muVec, lambda, muphisig$phiMat, muphisig$Sigma_Yi, verbose=FALSE
           )}, 
         y, MuPhiSig) 

  return(ret)
}


# Get the fixed components for a subject, with everything vectorized.
# The order of indices is (j, t) for mu, (t, j, k) for phi and (t, s, j, l) for Sigma_Y
MuPhiSigMult <- function(t, obsGrid, mu, phi, Sigma_Y) {

  #obsGrid <- signif(obsGrid, 14)
  K <- dim(phi)[3]
  phi <- aperm(phi, c(2, 1, 3))
  Sigma_Y <- aperm(Sigma_Y, c(3, 1, 4, 2)) # New Sigma_Y has entries (j, t, l, s)

  ret <- lapply(t, function(tvec) {
    #ind <- match(signif(tvec, 14), obsGrid)
    ind <- match(tvec, obsGrid)
    if (sum(is.na(ind)) != 0) {
      stop('Time point not found in obsGrid.')
    }

    muVec <- as.numeric(mu[, ind])
    phiMat <- matrix(phi[, ind, ], ncol=K)
    Sigma_Yi <- Sigma_Y[, ind, , ind, drop=FALSE]
    Sigma_Yi <- makeSquare(Sigma_Yi)
    
    return(list(muVec=muVec, 
                phiMat=phiMat, 
                Sigma_Yi=Sigma_Yi))
  })

  return(ret)
}


# GetIndCEScores <- function(yVec, muVec, lamVec, phiMat, Sigma_Yi, newyInd=NULL, verbose=FALSE) {

  # if (length(yVec) == 0) {
    # if (verbose) {
      # warning('Empty observation found, possibly due to truncation')
    # }
    # return(list(xiEst=matrix(NA, length(lamVec)), xiVar=matrix(NA, length(lamVec), length(lamVec)), fittedY=matrix(NA, 0, 0)))
  # }

# #
# ## When an individual has only one observation, the leave-one-out predicted Y is NA.
# #  if (length(yVec) == 1 && !is.null(newyInd)) {
# #    newPhi <- matrix(NA, ncol=length(lamVec))
# #    newMu <- NA
# #  }
# #
# #  if (!is.null(newyInd) && length(yVec) != 1) {
# #    # newy <- yVec[newyInd]
# #    newPhi <- phiMat[newyInd, , drop=FALSE]
# #    newMu <- muVec[newyInd]
# #
# #    yVec <- yVec[-newyInd]
# #    muVec <- muVec[-newyInd]
# #    phiMat <- phiMat[-newyInd, , drop=FALSE]
# #    Sigma_Yi <- Sigma_Yi[-newyInd, -newyInd, drop=FALSE]
# #  }
# #
# #  Lam <- diag(x=lamVec, nrow = length(lamVec))
# #  LamPhi <- Lam %*% t(phiMat)
# #  LamPhiSig <- LamPhi %*% solve(Sigma_Yi)
# #  xiEst <- LamPhiSig %*% matrix(yVec - muVec, ncol=1)
# #  xiVar <- Lam - LamPhi %*% t(LamPhiSig)
# #
# #
# #  fittedY <- if(is.null(newyInd)) 
# #    muVec + phiMat %*% xiEst else 
# #      newMu + newPhi %*% xiEst
# #
# #  ret <- list(xiEst=xiEst, xiVar=xiVar, fittedY=fittedY)
# #
# #  return(ret)

  # # Do all subscripting stuff in R
  # if (!is.null(newyInd)) {    
    # if (length(yVec) != 1){ 
      # newPhi <- phiMat[newyInd, , drop=FALSE]
      # newMu <- muVec[newyInd]
      # yVec <- yVec[-newyInd]
      # muVec <- muVec[-newyInd]
      # phiMat <- phiMat[-newyInd, , drop=FALSE]
      # Sigma_Yi <- Sigma_Yi[-newyInd, -newyInd, drop=FALSE]  
      # return ( GetIndCEScoresCPPnewInd( yVec, muVec, lamVec, phiMat, Sigma_Yi, newPhi, newMu) )
    # } else {   
      # # This should be an uncommon scenario
      # Lam <- diag(x=lamVec, nrow = length(lamVec))
      # LamPhi <- Lam %*% t(phiMat)
      # LamPhiSig <- LamPhi %*% solve(Sigma_Yi)
      # xiEst <- LamPhiSig %*% matrix(yVec - muVec, ncol=1)
      # xiVar <- Lam - LamPhi %*% t(LamPhiSig) 
      # return( list(xiEst=xiEst, xiVar = xiVar, fittedY=NA) )
    # }
  # } 
  # return( GetIndCEScoresCPP( yVec, muVec, lamVec, phiMat, Sigma_Yi) )
  # # Unfortunately function overloading is not yet available in Rcpp
  # # GetIndCEScoresCPPnewInd and GetIndCEScoresCPP are nearly identical.

# }

