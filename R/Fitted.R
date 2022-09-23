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

#' Obtains the fitted values from an `RFPCA` object
#' @param object An `RFPCA` object
#' @param K The number of components to apply. If missing, default to all components in `object`
#' @param grid The grid where the output function should be supported on
#' @param ... Not used
#' @export
#' @return A 3D array with entries corresponds to (i, j, t) supported on workGrid
fitted.RFPCA <- function(object, K, grid=c('work', 'obs'), ...) {

  grid <- match.arg(grid)
  if (missing(K)) {
    K <- object[['K']]
  }

  # buff <- 1e-15
  # ToutRange <- object[['optns']][['ToutRange']]
  # regGrid <- object[['regGrid']]
  # ind <- regGrid > ToutRange[1] - buff & regGrid < ToutRange[2] + buff
  # muRegTrunc <- object[['muReg']][, ind, drop=FALSE]
  if (grid == 'work') {
    mu <- object[['muWork']]
    phi <- object[['phi']]
    gridPt <- object[['workGrid']]
  } else if (grid == 'obs') {
    mu <- object[['muObsTrunc']]
    phi <- object[['phiObsTrunc']]
    gridPt <- object[['obsGridTrunc']]
  }

  res <- Fitted(object[['mfd']], mu, object[['xi']], phi, K)
  dimnames(res) <- list(i=rownames(object[['xi']]), 
                        j=rownames(mu), 
                        t=gridPt)
  res
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
