#' Riemannian Functional Principal Component Analysis
#'
#' FPCA for Riemannian manifold-valued functional data. The Riemannian (or multi-dimensional) functional data can be dense or sparse. 
#'
#' @param Ly A list of matrices, each being D by n_i containing the observed values for individual i. In each matrix, columes corresponds to different time points, and rows corresponds to different dimensions. 
#' @param Lt A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' 
#' @details Supported classes includes: 'Sphere' (default), 'Euclidean', 'SO', 'HS', 'L2', 'Dens'.
#'
#' Available control options are 
#' \describe{
#' \item{mfd}{A structure such as structure(1, class='CLASS'), where CLASS is one of the supported classes above. Takes precedence over mfdName. Default: structure(1, 'Sphere')}
#' \item{mfdName}{The name of a manifold. Supported values are 'Sphere', 'Euclidean', 'SO', 'HS', 'L2', 'Dens'. Default: 'Sphere'}
#' \item{dataType}{The type of design we have (usually distinguishing between sparse or dense functional data); 'Sparse', 'Dense', 'DenseWithMV', 'p>>n'. Default:  determine automatically based on 'IsRegular'}
#' \item{userBwMu}{The bandwidth for smoothing the mean function. Can be either a scalar specifying the bandwidth, or 'GCV' for generalized cross-validation. MUST BE SPECIFIED}
#' \item{userBwCov}{The bandwidth for smoothing the covariance function. Can be a scalar specifying the bandwidth. If userBwCov = 'GCV', then userBwCov will be set to twice the GCV-selected bandwidth for mu. MUST BE SPECIFIED}
#' \item{ToutRange}{Truncate the FPCA to be only within ToutRange. Default: c(-Inf, Inf)}
#' \item{npoly}{The degree of local polynomials for smoothing. Default: 1 (local linear)}
#' \item{nRegGrid}{The number of support points in each direction of covariance surface. Default: 51}
#' \item{kernel}{Smoothing kernel choice, common for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar". Default: "gauss"; dense data are assumed noise-less so no smoothing is performed. }
#' \item{error}{Assume measurement error in the dataset. The error is assumed to be isometric on the tangent space. Default: TRUE}
#' \item{maxK}{The maximum number of principal components to consider. Default: Inf if the smoothing method is used, and 30 for cross-sectional estimate}
#' \item{userSigma2}{The user-defined measurement error variance. A positive scalar. Default: `NULL`}
#' \item{methodMuCovEst}{The method to estimate the mean and covariance in the case of dense functional data; 'cross-sectional', 'smooth'. Default: 'cross-sectional'}
#' \item{methodXi}{The method to estimate the PC scores; 'CE' (Condit. Expectation), 'IN' (Numerical Integration). Default: 'CE' for sparse data and dense data with missing values, 'IN' for dense data.}
#' \item{obsGridOnly}{If TRUE, then assume the observation grids are regular, and use it as the regGrid/workGrid. This may speed up the grid convertion and eigendecomposition if length(obsGrid) is small. Default: TRUE if the Lt are regular and length(obsGrid) <= nRegGrid, and FALSE otherwise.}
#' }
#' @return A list containing the following fields:
#' \item{muReg}{A D by nRegGrid matrix containing the mean function estimate on regGrid.}
#' \item{muWork}{A D by nWorkGrid matrix containing the mean function estimate on workGrid.}
#' \item{muObs}{A D by nObsGrid matrix containing the mean function estimate on obsGrid.}
#' \item{cov}{An nWorkGrid by nWorkGrid by D by D array of the smoothed covariance surface.}
#' \item{covObs}{An nObsGrid by nObsGrid by D by D array of the smoothed covariance surface interpolated onto the obsGrid.}
#' \item{phi}{An nWorkGrid by D by \emph{K} array containing eigenfunctions supported on workGrid, where D is the ambient dimension.}
#' \item{phiObsTrunc}{A possibly truncated version of phi, supported on the truncated obsGrid}
#' \item{lambda}{A vector of length \emph{K} containing eigenvalues.}
#' \item{xi}{A \emph{n} by \emph{K} matrix containing the FPC estimates.} 
#' \item{sigma2}{Variance for measure error.}
#' \item{obsGrid}{The (sorted) grid points where all observation points are pooled.}
#' \item{regGrid}{A vector of length nRegGrid. The internal regular grid on which the eigen analysis is carried on.}
#' \item{workGrid}{Duplicates regGrid. A vector of length nWorkGrid. The internal regular grid on which the eigen analysis is carried on.}
#' \item{workGridTrunc}{A possibly truncated version of regGrid.}
#' \item{K}{Number of components returned}
#' \item{userBwMu}{The selected (or user specified) bandwidth for smoothing the mean function.}
#' \item{userBwCov}{The selected (or user specified) bandwidth for smoothing the covariance function.}
#' \item{mfd}{The manifold on which the analysis is performed.}
#' \item{optns}{A list of actually used options.}
#' @examples
#' # First simulate some data
#' set.seed(1)
#' n <- 50
#' m <- 20 # Number of different pooled time points
#' K <- 20
#' lambda <- 0.07 ^ (seq_len(K) / 2)
#' basisType <- 'legendre01'
#' sparsity <- 5:15
#' sigma2 <- 0.01
#' muList <- list(
#'   function(x) x * 2, 
#'   function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
#'   function(x) rep(0, length(x))
#' )
#' D <- length(muList)
#' pts <- seq(0, 1, length.out=m)
#' mfd <- structure(1, class='Euclidean')
#' mu <- Makemu(mfd, muList, c(rep(0, D - 1), 1), pts)
#'
#' # Generate noisy samples
#' samp <- MakeMfdProcess(mfd, n, mu, pts, K = K, lambda=lambda, basisType=basisType, sigma2=sigma2)
#' spSamp <- SparsifyM(samp$X, samp$T, sparsity)
#' yList <- spSamp$Ly
#' tList <- spSamp$Lt
#'
#' # Fit model
#' bw <- 0.2
#' kern <- 'epan'
#'
#' resEu <- RFPCA(yList, tList, 
#'   list(userBwMu=bw, 
#'        userBwCov=bw * 2, 
#'        kernel=kern, 
#'        maxK=K, 
#'        mfdName='euclidean', 
#'        error=TRUE))
#'
#' # Solid curve stands for the true mean and dashed for the estimated mean function.
#' matplot(pts, t(mu), type='l', lty=1)
#' matplot(pts, t(resEu$muObs), type='l', lty=2, add=TRUE)
#'
#' # Up to the 3rd principal components were well-estimated
#' plot(resEu$xi[, 3], samp$xi[, 3]) 
#' 
#' @export
RFPCA <- function(Ly, Lt, optns=list()) {

  optns <- SetOptionsRFPCA(Ly, Lt, optns)
  varnames <- rownames(Ly[[1]])

  # Evaluate all the option names
  for (optn in names(optns)) {
    assign(optn, optns[[optn]])
  }

  Ymat <- do.call(cbind, Ly)
  Tvec <- do.call(c, Lt)
  obsGrid <- sort(unique(Tvec))
  dimAmbient <- nrow(Ymat)
  dimTangent <- calcTanDim(mfd, dimAmbient=dimAmbient)
  ord <- order(Tvec)

  regGrid <- seq(min(obsGrid), max(obsGrid), length.out=nRegGrid)
  a <- ifelse(is.infinite(ToutRange[1]), min(obsGrid), ToutRange[1])
  b <- ifelse(is.infinite(ToutRange[2]), max(obsGrid), ToutRange[2])
  workGrid <- seq(a, b, length.out=nRegGrid)
  m <- length(workGrid)

  # Get bw
  if (dataType == 'Sparse') {
    if (userBwMu == 'GCV') {
      userBwMu <- GCVFrechetMeanCurve(mfd, Ly, Lt, kernel, 0, 'Sparse')$bOpt
      userBwCov <- 2 * userBwMu
    } 

    muReg <- frechetMeanCurve(mfd, userBwMu, kernel, xin=Tvec[ord], yin=Ymat[, ord],
                              xout=regGrid, npoly = npoly)
    muObs <- apply(muReg, 1, function(x) {
      ConvertSupport(as.numeric(regGrid), as.numeric(obsGrid), mu=x)
    })
    muObs <- apply(muObs, 1, project, mfd=mfd)
  } else if (dataType == 'Dense') {
    ## TODO: See which one is faster
    # muObs <- frechetMeanCurve(mfd, bw=(obsGrid[2] - obsGrid[1]) / 1e3, 
                              # 'epan', Tvec[ord], yin=Ymat[, ord], 
                              # xout=obsGrid, npoly = 0)
    muObs <- t(plyr::daply(data.frame(x=Tvec, t(Ymat)), 'x', function(dat) {
      y <- t(as.matrix(dat[-1]))
      c(frechetMean(mfd, y))
    }, .drop_o=FALSE))
    if (obsGridOnly) {
      muReg <- muObs
    } else {
      muReg <- t(apply(muObs, 1, function(x) {
        ConvertSupport(as.numeric(obsGrid), as.numeric(regGrid), mu=x)
      }))
      muReg <- apply(muReg, 2, project, mfd=mfd)
    }
  }

  if (isTRUE(all.equal(regGrid, workGrid))) {
    muWork <- muReg
  } else { # workGrid and regGrid are not the same
    muWork <- apply(muReg, 1, function(x) {
      ConvertSupport(as.numeric(regGrid), workGrid, mu=x)
    })
    muWork <- apply(muWork, 1, project, mfd=mfd)
  }

  rownames(muObs) <- varnames
  colnames(muObs) <- obsGrid
  rownames(muWork) <- varnames
  colnames(muWork) <- workGrid

  if (optns$meanOnly == TRUE) {
    return(list(muReg=muReg, regGrid=regGrid, 
                muWork=muWork, workGrid=workGrid))
  }

  if (meanOnly) {
    res <- list(muReg = muReg, muWork = muWork, muObs = muObs, 
                regGrid = regGrid, workGrid = workGrid, 
                obsGrid = obsGrid, 
                userBwMu=userBwMu, userBwCov=userBwCov, mfd = mfd, optns=optns)

    class(res) <- 'RFPCA'
    return(res)
  }

  # Log maps
  yListLog <- lapply(seq_along(Ly), function(i) {
    tt <- Lt[[i]]
    mu0 <- muObs[, match(tt, obsGrid), drop=FALSE]
    yy <- rieLog(mfd, mu0, Ly[[i]])
    yy
  })


  # Covariance
  if (dataType == 'Sparse') {
    covWork <- smoothCovM2(yListLog, Lt, matrix(0, nrow(muObs), ncol(muObs)), 
                           workGrid, userBwCov, kernel, error=error)
  } else if (dataType == 'Dense') {
    covObs <- csCovM(yListLog, Lt, error=FALSE)
    if (obsGridOnly) {
      covWork <- covObs
    } else {
      covWork <- array(
        apply(covObs, c(3, 4), 
              function(mat) {
                ConvertSupport(as.numeric(obsGrid), as.numeric(workGrid), Cov=mat, 
                               isCrossCov=TRUE)
              }), c(length(workGrid), length(workGrid), dimTangent, dimTangent)) 
    }
  }

  # Project to the correct tangent planes.
  covWork <- projectCov(mfd, covWork, muWork)


  # Truncate the analysis -- use only up to 
  # buff <- 1e-15
  # tIndTrunc <- regGrid >= ToutRange[1] - buff & regGrid <= ToutRange[2] + buff
  # workGrid <- regGrid[tIndTrunc]

  tmp <- lapply(seq_along(Lt), function(i) {
    tt <- Lt[[i]]
    yy <- yListLog[[i]]
    ind <- tt >= ToutRange[1] & tt <= ToutRange[2]
    list(t = tt[ind], y = yy[, ind, drop=FALSE])
  })

  LtTrunc <- sapply(tmp, `[[`, 't', simplify=FALSE)
  LyLogTrunc <- sapply(tmp, `[[`, 'y', simplify=FALSE)
  names(LtTrunc) <- names(LyLogTrunc) <- names(Lt)
  subInd <- sapply(LtTrunc, length) > 0
  if (!all(subInd)) {
    warning('Some subjects have no observations!')
    LtTrunc <- LtTrunc[subInd]
    LyLogTrunc <- LyLogTrunc[subInd]
  }

  TTruncInd <- obsGrid >= ToutRange[1] & obsGrid <= ToutRange[2]
  obsGridTrunc <- obsGrid[TTruncInd]
  mObsTrunc <- length(obsGridTrunc)

  muObsTrunc <- muObs[, TTruncInd, drop=FALSE]

  # Estimate the noise
  if (!is.null(userSigma2)) {
    sigma2 <- userSigma2
  } else if (error) {
    sigma2 <- EstSigma2(mfd, LyLogTrunc, LtTrunc, 
                        matrix(0, nrow(LyLogTrunc[[1]]), ncol(muObsTrunc)), 
                        covWork, workGrid, userBwCov, kernel, smooth=FALSE)
  } else {# error == FALSE
    sigma2 <- 0
  }

  # Eigenanalysis
  eig <- EigenAnalysis(covWork, workGrid, maxK, FVEthreshold, fastEig = fastEig, verbose=verbose)

  ## Normalization for functional observations
  functionalObs <- inherits(mfd, 'L2')
  if (functionalObs) {
    grids <- 1 / (dimTangent - 1) # grid size for s
    eig[['lambda']] <- eig[['lambda']] * grids
    eig[['phi']] <- eig[['phi']] / sqrt(grids)
  }
  lam <- eig[['lambda']]
  K <- length(lam)
  phi <- array(eig[['phi']], c(m, dimTangent, K))
  covFitted <- eig[['fittedCov']]

  # Scores
  if (methodMuCovEst == 'smooth') {
    covObs <- apply(eig[['fittedCov']], c(3, 4), 
                    function(mat) {
                      ConvertSupport(as.numeric(workGrid), as.numeric(obsGridTrunc), Cov=mat, 
                                     isCrossCov=TRUE)
                    }) 
    covObs <- array(covObs, c(mObsTrunc, mObsTrunc, dimTangent, dimTangent))
    covObs <- projectCov(mfd, covObs, muObsTrunc)
  }

  if (obsGridOnly) {
    phiObsTrunc <- phi
  } else {
    phiObsTrunc <- apply(phi, c(2, 3), function(phi) {
                      ConvertSupport(as.numeric(workGrid), as.numeric(obsGridTrunc), mu=phi)
                    })
  }
  phiObsTrunc <- projectPhi(mfd, phiObsTrunc, muObsTrunc)

  if (methodXi == 'CE') {
    CE <- CEScores(LyLogTrunc, 
                   LtTrunc, 
                   list(), 
                   matrix(0, dimTangent, mObsTrunc), 
                   obsGridTrunc, 
                   covObs, lam, phiObsTrunc, sigma2)
    xi <- t(do.call(cbind, CE['xiEst', ]))
    if (functionalObs) {
      stop('Not implemented yet')
    }
  } else if (methodXi == 'IN') {
    ## The time and s points must be regular for the time being
    ylogTrunc <- aperm(simplify2array(LyLogTrunc), c(3, 2, 1)) # Indices are (i, t, j)
    dims <- dim(ylogTrunc)
    n <- dims[1]
    yMat2 <- matrix(ylogTrunc, n, dims[2] * dims[3])
    phiMat <- matrix(phiObsTrunc, ncol=dim(phiObsTrunc)[3]) # Row index are (t, j)
    if (length(obsGridTrunc) == 1) {
      gridt <- 1
    } else {
      gridt <- mean(diff(obsGridTrunc))
    }
    xi <- yMat2 %*% phiMat * gridt
    if (functionalObs) {
      xi <- xi * grids
    }
  }

  rownames(xi) <- names(LtTrunc)
  colnames(xi) <- paste0('xi', seq_len(ncol(xi)))

  res <- list(muReg = muReg, muWork = muWork, 
              muObs = muObs, muObsTrunc = muObsTrunc, 
              cov = covWork, # On workGrid
              covFitted = covFitted, # On workGrid
              phi = phi, # On workGrid
              covObs = covObs, 
              phiObsTrunc = phiObsTrunc, 
              lam = lam, xi=xi,
              sigma2 = sigma2, 
              regGrid = regGrid, workGrid = workGrid, 
              obsGrid = obsGrid, obsGridTrunc=obsGridTrunc,
              K = K, userBwMu=userBwMu, userBwCov=userBwCov, mfd = mfd, optns=optns)

  class(res) <- 'RFPCA'

  res
}


# flip the signs of components
flipComponents <- function(res, flipInd) {

  if (missing(flipInd) || length(flipInd) == 0) {
    stop('flipInd has something wrong')
  }

  res$phi[, , flipInd] <- -res$phi[, , flipInd]
  res$phiObs[, , flipInd] <- -res$phiObs[, , flipInd]
  res$xi[, flipInd] <- -res$xi[, flipInd]

  res
}

