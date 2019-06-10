# Multivariate entrywise 2D smoothing. This is equivalent to the global smoothing using matrix Frobenius norm. Note that if the input y are not on the same tangent plane, the output covariance may not on a tangent plane.
# mu: Supported on the unique time points
# error: Whether measurement errors are assumed. 
smoothCovM2 <- function(yList, tList, mu, regGrid, bwCov, kernel_type, error=TRUE) {

  n <- length(yList)
  nT <- ncol(mu) # number of time points
  d <- nrow(yList[[1]])
  m <- length(regGrid)

  multyList <- lapply(seq_len(d), function(j) lapply(yList, `[`, j, ))

  covR <- 
    sapply(seq_len(d), function(j2) {
    sapply(seq_len(d), function(j1) {
      rcov <- fdapace:::GetRawCrCovFuncFunc(
        multyList[[j1]], tList, mu[j1, ], 
        multyList[[j2]], tList, mu[j2, ]
        )
      if (error) {
        ind <- rcov$tpairn[, 1] != rcov$tpairn[, 2]
      } else {
        ind <- rep(TRUE, nrow(rcov$tpairn))
      }
      # browser()
      res <- matrix(Lwls2D(bwCov, kernel_type, 
                           rcov$tpairn[ind, , drop=FALSE], 
                           rcov$rawCCov[ind], 
                           xout1=regGrid, xout2=regGrid, crosscov=TRUE), m, m)
      res
    }, simplify='array')
    }, simplify='array')

  covR
}

# Scalar version: Project the covariance surface to the tangent spaces at mu
projectCovS <- function(mfd=structure(1, class='Sphere'), covR, mu) {

  stopifnot(dim(covR)[1] == dim(covR)[2] && dim(covR)[1] == length(mu))
  mu <- matrix(mu)

  projM <- projectTangent(mfd, mu, projMatOnly=TRUE)
  covR <- projM %*% covR %*% t(projM)

  covR
}


# Longitudinal version: Project the covariance surface to the tangent spaces at mu
projectCov <- function(mfd=structure(1, class='Sphere'), covR, mu) {

  stopifnot(dim(covR)[1] == dim(covR)[2] && dim(covR)[1] == ncol(mu))
  m <- dim(covR)[1]

  projM <- lapply(seq_len(ncol(mu)), 
                  function(i) projectTangent(mfd, mu[, i], projMatOnly=TRUE))

  for (j1 in seq_len(m)) {
  for (j2 in seq_len(m)) {
    covR[j1, j2, , ] <- as.matrix(projM[[j1]] %*% covR[j1, j2, , ] %*% Matrix::t(projM[[j2]]))
  }
  }

  covR
}


# Project the eigenfunctions to the tangent spaces at mu
# Entries of phi are [t, j, k]
projectPhi <- function(mfd=structure(1, class='Sphere'), phi, mu) {

  dimTangent <- dim(phi)[2]
  dimAmbient <- nrow(mu)
  m <- dim(phi)[1]
  
  stopifnot(dimTangent == calcTanDim(mfd, dimAmbient=dimAmbient) && 
            dim(phi)[1] == ncol(mu))

  projM <- lapply(seq_len(ncol(mu)), 
                  function(i) projectTangent(mfd, mu[, i], projMatOnly=TRUE))

  for (j in seq_len(m)) {
    phi[j, , ] <- as.matrix(projM[[j]] %*% phi[j, , ])
  }

  phi
}


# Cross-sectional covariance estimate on a manifold. Assume the elements in yList are all centered.
csCovM <- function(yList, tList, error=TRUE) {

  if (!error) {
    # browser()
    yAr <- simplify2array(yList)
    dims <- dim(yAr)
    n <- dims[3]
    yMat <- matrix(yAr, dims[1] * dims[2], dims[3])
    res <- array(tcrossprod(yMat) / (n - 1), c(dims[1], dims[2], dims[1], dims[2]))
    res <- aperm(res, c(2, 4, 1, 3))
  } else if (error) {
    stop('Cross-sectional covariance estimate with error is not implemented')
  }

  res
}

# # TODO: get the local windows, or indices in the local windows
# # 
# # Obtain smooth estimate of the covariance function on tangent spaces.
# # FOR NOW, assume mu is on the regGrid
# # yList: A list of function values in matrices, whose columns correspond to different measurements.
# # tList: A list of time points in vectors
# # mu: A matrix of mean curve
# # regGrid: Output support of of the covariance function
# # bwCov: A scalar bandwidth for smoothing the covariance
# # kernel_type: Kernel for smoothing, only Epanechnikov and Gaussian kernels are supported.
# # error: Whether to remove the diagonal observation time points. Same as whether measurement error is assumed.

# smoothCovM <- function(yList, tList, mu, regGrid, bwCov, kernel_type, error=TRUE) {

  # n <- length(yList)
  # nT <- ncol(mu) # number of time points
  # p <- nrow(mu)
  # stopifnot(nT == length(regGrid))
  # m <- length(regGrid)
  # eps <- 1e-10

  # # Interpolate the mean onto outGrid. SKIP NOW

  # # Get all log data and raw covariances
  # ally <- do.call(cbind, yList)
  # allt <- do.call(c, tList)
  # allID <- rep(seq_len(n), sapply(tList, length))

  # # For each regGrid, calculate the logarithm map if the observation is within a local window
  # if (kernel_type == 'gauss') {
    # bwWindow <- Inf
  # } else {
    # bwWindow <- bwCov
  # }

  # # For matrix logDat: t0 is the time for the base point, t is the time for the observation
  # logDat <- do.call(rbind, 
    # sapply(regGrid, function(tt) {
             
             # # browser()
      # # Get the time window
      # ind <- abs(tt - allt) <= eps + bwWindow
      # subt <- allt[ind]
      # suby <- ally[, ind, drop=FALSE]
      # tind <- match(tt, regGrid)

      # mfd <- structure(1, class='Sphere2D')
      # suby <- rieLog(mfd, mu[, tind, drop=FALSE], suby) # take log here
      # subID <- allID[ind]
      # cbind(y = t(suby), ID = subID, t = subt, t0 = tt)

    # }, simplify=FALSE)
  # )

  # # An array to store the results. Indices are (t, s, j1, j2), where (j1, j2) are entries in the pointwise covariance matrix.
  # res <- array(as.numeric(NA), dim=c(m, m, p, p))

  # for (indt0 in seq_along(regGrid)) {
  # for (inds0 in seq_along(regGrid)) {
    # if (indt0 > inds0) { # fill by symmetry
      # res[indt0, inds0, , ] <- t(res[inds0, indt0, , ])
    # } else { # calculate
      # t0 <- regGrid[indt0]
      # s0 <- regGrid[inds0]
      # tSamp <- logDat[logDat[, 't0'] == t0, , drop=FALSE]
      # sSamp <- logDat[logDat[, 't0'] == s0, , drop=FALSE]
      # useID <- intersect(unique(tSamp[, 'ID']), unique(sSamp[, 'ID']))

      # tmp <- sapply(useID, function(ID) {
        # ttt <- tSamp[tSamp[, 'ID'] == ID, , drop=FALSE]
        # sss <- sSamp[sSamp[, 'ID'] == ID, , drop=FALSE]
        # ycols <- seq_len(which(colnames(ttt) == 'ID') - 1)
        # tPairs <- as.matrix(expand.grid(t=ttt[, 't'] - t0, s=sss[, 't'] - s0))

        # # Get the raw covariances as matrices
        # tmp <- outer(ttt[, ycols, drop=FALSE], sss[, ycols, drop=FALSE])
        # rcov <- matrix(aperm(tmp, c(1, 3, 2, 4)), nrow = nrow(ttt) * nrow(sss))
        # list(yy = rcov, tPairs = tPairs)
      # }, simplify=FALSE)
      # rcov <- list(
        # yy = do.call(rbind, lapply(tmp, `[[`, 'yy')), 
        # tPairs = do.call(rbind, lapply(tmp, `[[`, 'tPairs'))
        # ) 

      # # Get weights based on tPairs
      # wt <- kernelWeight2D(rcov[['tPairs']], bwCov, kernel_type, eps)

      # # Smooth with each window
      # # browser()
      # tsMat <- smoothRawCov1(rcov[['tPairs']], rcov[['yy']], wt)

      # res[indt0, inds0, , ] <- tsMat
    # }
  # }
  # }

  # res

# }

# Get kernel weights

kernelWeight2D <- function(xin, h, kernel_type, teval=c(0, 0), eps=1e-10) {

  xin <- matrix(xin, ncol=2)
  N <- nrow(xin)
  xin <- xin - matrix(teval, nrow=N, ncol=2, byrow=TRUE)

  if (substr(kernel_type, 1, 4) != 'gaus') {
    stopifnot(all(abs(xin) <= h + 2 * eps))
  }

  # The normalizing constants does not matter here. Implemented just Epanechnikov and Gaussian kernel now.
  weight <- switch(
    kernel_type, 
    'epan' = (1 - (xin[, 1] / h)^2) * (1 - (xin[, 2] / h)^2),
    # 'rect' = ,
    'gauss' = exp(- (xin[, 1] / h)^2 / 2 - (xin[, 2] / h)^2 / 2)#,
    # 'gausVar' = ,
    # 'quar' = ,
  )

  weight / sum(weight)
}

# combLogDat <- function(logDatList) {

# }


# Get the smoothed covariance matrix at a certain time point based on a local linear smoother, given raw data mapped onto the tangent planes and kernel weights.
# The entrywise method is mathematically equivalent to the global method. Use the former.
# xin: n by 2 matrix
# yin: n by p^2 matrix
# win: a vector of length n
smoothRawCov1 <- function(xin, yin, win) {
  X <- cbind(1, xin)
  W <- diag(win, length(win))
  # wt <- solve(crossprod(X, W %*% X), crossprod(X, W))[1, ]
  wt <- (MASS::ginv(crossprod(X, W %*% X)) %*% crossprod(X, W))[1, ]

  linComb(wt, yin, makeMat=TRUE)
}


# a linear combination of matrices
# a: Coefficients in the linear combination
# M: A matrix with n rows, each representing a matrix
# makeMat: Whether to make the output a square matrix.
linComb <- function(a, M, makeMat=FALSE) {
  stopifnot(is.vector(a))
  res <- as.numeric(a %*% M)
  if (makeMat) {
    res <- makeSquare(res)
  }

  res
}

makeSquare <- function(v) {
  l <- length(v)
  p <- floor(sqrt(l))
  stopifnot(p == sqrt(l))
  res <- matrix(v, p, p)

  res
}


