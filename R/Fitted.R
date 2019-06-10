# Return: A 3D array with entries corresponds to (i, j, t)
Fitted <- function(mfd=structure(1, class='Sphere'), muM, xi, phiM, K) {

  if (missing(K)) {
    K <- ncol(xi)
  }

  n <- nrow(xi)
  m <- dim(phiM)[1] # number of time points
  p <- dim(phiM)[2] # ambient dimension

  xi <- xi[, seq_len(K), drop=FALSE]
  phiM <- phiM[, , seq_len(K), drop=FALSE]
  phiMat <- matrix(phiM, m * p, K)

  xiphi <- tcrossprod(xi, phiMat)

  res <- sapply(seq_len(n), function(i) {
    tmp <- sapply(seq_len(m), function(j) {
      rieExp(mfd, 
             muM[, j, drop=FALSE], 
             t(matrix(xiphi[i, ], m, p))[, j, drop=FALSE])
    })
  }, simplify='array')
  
  aperm(res, c(3, 1, 2))
}


#' Get the fitted trajectories for an RFPCA object
#'
#' @param K The number of components to use
#' @param grid Either `obs` for obsGrid or `work` for workGrid
#' @export
fitted.RFPCA <- function(object, K, grid=c('obs', 'work'), ...) {

  grid <- match.arg(grid)
  res <- object
  if (missing(K)) {
    K <- res[['K']]
  }

  # buff <- 1e-15
  # ToutRange <- res[['optns']][['ToutRange']]
  # regGrid <- res[['regGrid']]
  # ind <- regGrid > ToutRange[1] - buff & regGrid < ToutRange[2] + buff
  # muRegTrunc <- res[['muReg']][, ind, drop=FALSE]

  if (grid == 'obs') {
    mu <- res[['muObs']]
    phi <- res[['phiObsTrunc']]
  } else if (grid == 'work') {
    mu <- res[['muWork']]
    phi <- res[['phi']]
  }

  Fitted(res[['mfd']], mu, res[['xi']], phi, K)
}


PositiveOrth <- function(X) {

  X[X < 0] <- 0
  X[X > 1] <- 1

  dd <- dim(X)
  if (length(dd) == 2) { # entries corresponds to (j, t)
    X <- apply(X, 2, Normalize)
    res <- array(X, dd)
  } else if (length(dd) == 3) { # entries corresponds to (i, j, t)
    X <- apply(X, c(1, 3), Normalize)
    res <- aperm(X, c(2, 1, 3))
    # res <- X
  }

  res
}
