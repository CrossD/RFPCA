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



#' @export
# If type is 'traj', then predict the trajectories on the obsGrid using the first K components
predict.RFPCA.HS <- function(object, newLy, newLt, sigma2, K, xiMethod, type=c('xi', 'traj'), ...) {

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

#' Obtains the predicted value for a set of new sparse observations
#' @param object An `RFPCA` object
#' @param newLy A list containing the new response values; each entry corresponds to a subject and should be a matrix in which each column is an observation. 
#' @param newLt A list containing the new time points; each entry corresponds to a subject and should be a vector.
#' @param sigma2 The noise variance to use for the prediction. Serves as a regularization parameter. If left as `NULL`, use the one chosen when fitting `object`.
#' @param K The number of components to apply. If missing, default to all components in `object`
#' @param xiMethod Either to use the conditional expectation (`'CE'`) or the integration method (`'IN'`) to predict
#' @param type Either to predict the FPC scores (`'xi'`) or the underlying trajectories (`'traj'`)
#' @export
#' @return A 2D array of predicted scores with entries corresponds to (i, k)
predict.RFPCA <- function(object, newLy, newLt, 
                          sigma2 = NULL, K, xiMethod = c('CE', 'IN'), 
                          type=c('xi', 'traj'), ...) {

  xiMethod <- match.arg(xiMethod)

  if (missing(K)) {
    K <- object[['K']]
  }

  if(! all.equal(sapply(newLt, length), sapply(newLy, ncol) ) ){
    stop('The size of the vectors in newLt and newLy differ. They must be equal.')
  } 
  
  # TODO: check missing values  

  if(is.null(sigma2)) {
    sigma2 <- ifelse( !is.null(object$rho), object$rho, 
                      ifelse( !is.null(object$sigma2), object$sigma2,
                              ifelse( !object$optns$error, 0, stop('sigma2 cannot be determined.'))))  
  } else {
    if(!is.numeric(sigma2) ) { 
      stop('sigma2 is not numeric, we sure about this? :D')
    }
  }
  
  if(K > object$K) {
    stop( paste0( collapse = '', 'You cannot get FPC scores for more components than what is already available. (', object$selectK ,').' ))
  } 
  
  if( (object$optns$dataType == 'Sparse') && (xiMethod == 'IN')){
    stop( 'Trapezoid Numerical intergration (IN) is invalid for sparse data.')
  }
  
  ToutRange <- range(object$workGrid)
  mfd <- object$optns$mfd
  dimAmbient <- nrow(newLy[[1]])
  dimTangent <- calcTanDim(mfd, dimAmbient=dimAmbient)
  funcObs <- inherits(mfd, 'L2')

  muWork <- object$muWork
  workGrid <- object$workGrid
  obsGrid <- object$obsGrid

  obsGridNew <- sort(unique(unlist(newLt)))
  obsGridNew <- obsGridNew[obsGridNew >= min(ToutRange) & 
                           obsGridNew <= max(ToutRange)] # Observations out of the training time domain are not used.

                         # browser()
  useGrid <- 
    if (length(obsGridNew) == length(obsGrid) && 
        max(abs(obsGridNew - obsGrid)) <= 1e-14) {
      'obsGrid'
    } else if (length(obsGridNew) == length(workGrid) && 
               max(abs(obsGridNew - workGrid)) <= 1e-14) {
      'workGrid'
    } else {
      'convert'
    }

  if (useGrid == 'obsGrid') {
    muObsNew <- object$muObsTrunc
  } else if (useGrid == 'workGrid') {
    muObsNew <- object$muWork
  } else if (useGrid == 'convert') {
    muObsNew <- apply(muWork, 1, function(x) {
      ConvertSupport(workGrid, obsGridNew, mu=x)
    })
    muObsNew <- apply(muObsNew, 1, project, mfd=mfd)
  }


  # perform truncation
  tmp <- lapply(seq_along(newLt), function(i) {
    tt <- newLt[[i]]
    ind <- tt >= ToutRange[1] & tt <= ToutRange[2]
    tt <- tt[ind]
    yy <- newLy[[i]][, ind, drop=FALSE]

    mu0 <- muObsNew[, match(tt, obsGridNew), drop=FALSE]
    yy <- rieLog(mfd, mu0, yy)
    list(t = tt, y = yy)
  })

  LtTrunc <- sapply(tmp, `[[`, 't', simplify=FALSE)
  LyLogTrunc <- sapply(tmp, `[[`, 'y', simplify=FALSE)
  names(LtTrunc) <- names(LyLogTrunc) <- names(newLt)

  subInd <- sapply(LtTrunc, length) > 0
  if (!all(subInd)) {
    warning('Some subjects have no observations!')
    LtTrunc <- LtTrunc[subInd]
    LyLogTrunc <- LyLogTrunc[subInd]
  }

  mObsTrunc <- length(obsGridNew)

  # browser()

  if (useGrid == 'obsGrid') {
    covUse <- object[['covObs']]
    phiUse <- object[['phiObsTrunc']]
  } else if (useGrid %in% c('workGrid', 'convert')) {
    covUse <- object[['covFitted']]
    phiUse <- object[['phi']]
  }

  covObs <- apply(covUse, 
                  c(3, 4), 
                  function(mat) {
                    if (useGrid == 'convert') {
                      ConvertSupport(workGrid, obsGridNew, Cov=mat, 
                                     isCrossCov=TRUE)
                    } else {
                      mat
                    }
                  })
  covObs <- array(covObs, c(mObsTrunc, mObsTrunc, dimTangent, dimTangent))
  covObs <- projectCov(mfd, covObs, muObsNew)

  phiObsNew <- apply(phiUse[, , seq_len(K), drop=FALSE], c(2, 3), function(phi) {
                    if (useGrid == 'convert') {
                      ConvertSupport(workGrid, obsGridNew, mu=phi)
                    } else {
                      phi
                    }
                  })
  phiObsNew <- array(phiObsNew, c(mObsTrunc, dimTangent, K))
  phiObsNew <- projectPhi(mfd, phiObsNew, muObsNew)

  if (xiMethod == 'CE') {
    if (type == 'traj') stop('type=traj is to be implemented for xiMethod=CE')

    CE <- CEScores(LyLogTrunc, 
                   LtTrunc, 
                   list(), 
                   matrix(0, dimTangent, mObsTrunc), 
                   obsGridNew, 
                   covObs, object[['lam']], phiObsNew, sigma2)

    xi <- t(do.call(cbind, CE['xiEst', ]))
    rownames(xi) <- names(LtTrunc)
    colnames(xi) <- paste0('xi', seq_len(ncol(xi)))

    res <- xi

  } else if (xiMethod == 'IN') {
    ymat <- aperm(simplify2array(LyLogTrunc), c(3, 1, 2))
    res <- PredictIn(mfd, ymat, obsGridNew, muObsNew, phiObsNew, K, type, funcObs)
  }
  res
}
