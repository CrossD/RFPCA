# muM The mean function on obsGrid
# newY A 3D array with entries corresponds to (i, j, t), for t on obsGrid
# Return: A 3D array with entries corresponds to (i, j, t)
# For now, assume obsGrid is equally spaced.
PredictIn <- function(mfd=structure(1, class='Sphere'), newY, obsGrid, muM, phiM, K, type, funcObs) {

  n <- dim(newY)[1]
  m <- dim(phiM)[1] # number of time points
  p <- dim(phiM)[2] # ambient dimension
  J <- dim(phiM)[3]

  if (missing(K)) {
    K <- J
  }

  # Project X onto phi
  phiM <- phiM[, , seq_len(K), drop=FALSE]
  phiMat <- matrix(phiM, m * p, K)
  Ymat <- matrix(aperm(newY, c(1, 3, 2)), n, m * p)
  if (m == 1) {
    gridSize <- 1
  } else {
    gridSize <- mean(diff(obsGrid))
  }
  xi <- Ymat %*% phiMat * gridSize
  if (funcObs) {
    xi <- xi / (p - 1)
  }

  if (type == 'xi') {
    res <- xi
  } else if (type == 'traj') {
    xiphi <- tcrossprod(xi, phiMat)
    res <- sapply(seq_len(n), function(i) {
      tmp <- sapply(seq_len(m), function(j) {
        V <- t(matrix(xiphi[i, ], m, p))
        rieExp(mfd, 
               muM[, j, drop=FALSE], 
               V[, j, drop=FALSE])
      })
    }, simplify='array')
    res <- aperm(res, c(3, 1, 2))
  }
  res
}


# If type is 'traj', then predict the trajectories on the obsGrid using the first K components
predict.RFPCA <- function(object, newLy, newLt, sigma2, K, xiMethod, type=c('xi', 'traj'), ...) {

  if (!missing(sigma2)) {
    stop('sigma2 not implemented')
  }
  type <- match.arg(type)
  mfd <- object[['mfd']]
  funcObs <- inherits(mfd, 'L2')

  xiMethod <- object[['optns']][['methodXi']]
  muObs <- object[['muObs']]
  phiObs <- object[['phiObsTrunc']]
  obsGrid <- object[['obsGrid']]

  if (missing(K)) {
    K <- object[['K']]
  }

  # Take log maps
  yListLog <- lapply(seq_along(newLy), function(i) {
    tt <- newLt[[i]]
    mu0 <- muObs[, match(tt, obsGrid), drop=FALSE]
    yy <- rieLog(mfd, mu0, newLy[[i]])
    yy
  })

  if (xiMethod == 'IN') {
    ymat <- aperm(simplify2array(yListLog), c(3, 1, 2))
    res <- PredictIn(mfd, ymat, obsGrid, muObs, phiObs, K, type, funcObs)
  } else if (xiMethod == 'CE') {
    stop('Not implemented yet')
    # res <- PredictCE(mfd, ymat, muObs, phiObs, K)
  }

  res
}
