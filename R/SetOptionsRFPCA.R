# mfdName 'Sphere', 'SO', 'Euclidean', etc
# mfd An alternative way to specify manifold. Takes precedence over mfdName
# obsGridOnly If TRUE, then assume the observation grids are regular, and use it as the regGrid/workGrid. This may speed up the grid convertion and eigendecomposition if length(obsGrid) is small. Default: Check Lt: TRUE if all Lt's are regular and length(obsGrid) <= nRegGrid, and FALSE otherwise.
# meanOnly If TRUE, then just estimate the mean function. Default: FALSE
SetOptionsRFPCA <- function(Ly, Lt, optns) {

  # Note: Must modify this line when adding a new option.
  allOptnNames <- c('mfdName', 'mfd', 'userBwMu', 'userBwCov', 'ToutRange', 'npoly', 
                    'nRegGrid', 'kernel', 'error', 'maxK', 'userSigma2', 'dataType', 
                    'methodMuCovEst', 'methodXi', 'obsGridOnly', 'meanOnly', 
                    'fastEig', 'FVEthreshold', 'verbose')

  invalidNames <- setdiff(names(optns), allOptnNames)

  if (any(nchar(names(optns)) == 0)) {
    stop('All specified options in `optns` must be named')
  }

  if (length(invalidNames) > 0) {
    stop(sprintf('Invalid option name(s): %s', paste(invalidNames)))
  }

  # Load all optns
  for (optn in allOptnNames) {
    expr <- parse(text=sprintf("%s <- optns[['%s']]", optn, optn))
    eval(expr)
  }

  if (is.null(dataType)) { #do we have dataType or sparse functional data
    dataType <- fdapace:::IsRegular(Lt)
  }


  if (is.null(methodMuCovEst)) {
    if (dataType == 'Sparse'){
      methodMuCovEst <- 'smooth'
      methodXi <- 'CE'
    } else {
      methodMuCovEst <- 'cross-sectional'
      methodXi <- 'IN'
    }
  }


  # Set defaults for options
  if (is.null(mfd)) {
    if (is.null(mfdName)) {
      mfdName <- 'Sphere'
    } 
    mfdName <- tolower(mfdName)
    if (mfdName == 'sphere') {
      mfd <- structure(1, class='Sphere')
    } else if (mfdName == 'euclidean') {
      mfd <- structure(1, class='Euclidean') 
    } else if (mfdName == 'so') {
      mfd <- structure(1, class='SO')
    } else if (mfdName == 'logeu') {
      mfd <- structure(1, class=c('LogEu', 'SPD'))
    } else if (mfdName == 'affsym') {
      mfd <- structure(1, class=c('AffSym', 'SPD'))
    }
  }
  mfdName <- tolower(class(mfd)[1]) # In case mfd is specified but not mfdName

  if (methodMuCovEst == 'smooth') {
    if (is.null(userBwMu)) {
      stop('Specify bandwidth for smoothing mu')
    }

    if (is.null(userBwCov)) {
      stop('Specify bandwidth for smoothing cov')
    }
  }

  if (is.null(ToutRange)) {
    ToutRange <- c(-Inf, Inf)
  }

  if (is.null(npoly)) {
    npoly <- 1
  }

  if (is.null(nRegGrid)) {
    nRegGrid <- 51
  }

  if (is.null(kernel)) {
    kernel <- 'epan'
  }

  if (is.null(error)) {
    if (dataType == 'Sparse') {
      error <- TRUE
    } else if (dataType == 'Dense') {
      error <- FALSE
    }
  }

  if (is.null(maxK)) {
    if (methodMuCovEst == 'smooth') {
      maxK <- Inf
    } else if (methodMuCovEst == 'cross-sectional') {
      maxK <- 30
    }
  }

  if (is.null(userSigma2)) {
    userSigma2 <- NULL
  }

  if (is.null(obsGridOnly)) {
    if (dataType == 'Dense' && 
        all(sapply(Lt, function(tt) all(diff(tt) == tt[2] - tt[1]))) &&
        length(Lt[[1]]) <= nRegGrid) {
      obsGridOnly <- TRUE
    } else {
      obsGridOnly <- FALSE
    }
  }
  
  if (obsGridOnly) {
    nRegGrid <- length(Lt[[1]])
  }

  if (is.null(meanOnly)) {
    meanOnly <- FALSE
  }

  if (is.null(FVEthreshold)) {
     FVEthreshold <- 1
  }

  if (is.null(verbose)) {
    verbose <- TRUE
  }

  if (is.null(fastEig)) {
    fastEig <- dataType == 'Dense'
  }

  # Get the returned list
  retExpr <- parse(
    text= sprintf('list(%s)', 
                  paste(sprintf('%s = %s', allOptnNames, allOptnNames), 
                                 collapse = ',\n'))
  )

  ret <- eval(retExpr)
  ret
}


AppendOptionsRFPCA <- function(optns, ...) {
  newOptns <- list(...)
  # stopifnot(all(names(newOptns) %in% names(optns)))
  optns[names(optns) %in% names(newOptns)] <- NULL
  optns <- c(optns, newOptns)
  optns
}
