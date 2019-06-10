#' Create a mean function 
#'
#' Generate a mean function through the exponential map. Specify a basepoint \eqn{p_0} and \eqn{V_\mu(t)}. Then \eqn{\mu=\exp_{p_0}V_{\mu(t)}
#' @param mfd An object whose class specifies the manifold to use
#' @param VtList A list of functions. The jth function takes a vector of time, and return the corresponding V(t) in the jth entry
#' @param p0 The basepoint, which will be converted to a vector
#' @param pts The time points to generate the mean function
#' @export
Makemu <- function(mfd, VtList, p0, pts=seq(0, 1, length.out=50)) {

  # UseMethod('Makemu', mfd)
  p0 <- c(p0)
  Vmat <- matrix(sapply(VtList, function(f) f(pts)), ncol=length(VtList))
  res <- rieExp(mfd, p0, t(Vmat))

  res
}

Normalize <- function(v, tol=1e-10) {
  v <- as.numeric(v)
  d <- length(v)
  vNorm <- as.numeric(sqrt(crossprod(v)))
  res <- if (!is.na(vNorm) && vNorm <= tol) {
    rep(0, d)
  } else {
    v / vNorm
  }

  res
}


## Project p1 onto p2. If rej is TRUE return the rejection instead
Proj <- function(p1, p2, rej=FALSE, tol=1e-10) {
  p1 <- as.numeric(p1)
  p2 <- Normalize(as.numeric(p2), tol)

  if (sum(p2^2) == 0) {
    stop('p2 cannot be 0')
  }

  proj <- c(crossprod(p1, p2)) * p2
  res <- if (!rej) {
    proj
  } else { # Return the rejection
    p1 - proj
  }

  as.numeric(res)
}


## Generate mean function by exponential map. 
# VtList: a list of coordinate functions for V(t). 
# Makemu.Sphere2D <- function(mfd, VtList, p0, pts=seq(0, 1, length.out=50)) {

  # Vmat <- sapply(VtList, function(f) f(pts))
  # res <- rieExp(mfd, p0, t(Vmat))

  # res
# }


# Makemu.Sphere <- function(mfd, VtList, p0, pts=seq(0, 1, length.out=50)) {

  # Vmat <- sapply(VtList, function(f) f(pts))
  # res <- rieExp(mfd, p0, t(Vmat))

  # res
# }


MakeRotMat <- function(p1, p2, tol=1e-10) {
  d <- length(p1)
  if (length(p1) != length(p2)) {
    stop('dimensions of p1 and p2 must agree')
  }
  if (sum(abs(p1 - Normalize(p1))) > tol || 
      sum(abs(p2 - Normalize(p2))) > tol) {
    stop('p1 and p2 needs to be have norm 1')
  }

  Sphere <- structure(1, class='Sphere')
  psi <- as.numeric(distance(Sphere, matrix(p1), matrix(p2)))
  p1t <- Normalize(Proj(p1, p2, rej=TRUE))
  R <- diag(1, nrow=d) + 
    sin(psi) * (tcrossprod(p2, p1t) - tcrossprod(p1t, p2)) + 
    (cos(psi) - 1) * (tcrossprod(p2, p2) + tcrossprod(p1t, p1t))

  R
}


## Generate eigenfunctions V(t) from a certain basis around a mean curve.
# mu needs to correspond to pts.
MakephiS <- function(K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
                     type=c('cos', 'sin', 'fourier', 'legendre01')) {

  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }
  p <- nrow(mu)
  if (p <= 1) {
    stop('mu must have more than 1 rows')
  }
  m <- length(pts)
  if (ncol(mu) != m) {
    stop('mu must have the same number of columns as the length of pts')
  }

  if (min(pts) < 0 || max(pts) > 1) {
    stop('The range of pts should be within [0, 1]')
  }

  ptsSqueeze <- do.call(c, lapply(seq_len(p - 1), function(i) {
                          pts + (max(pts) + mean(diff(pts))) * (i - 1)
                        }))
  ptsSqueeze <- ptsSqueeze / max(ptsSqueeze)
  
  phi1 <- rbind(fdapace:::CreateBasis(K, ptsSqueeze, type),
                matrix(0, m, K)) / sqrt(p - 1) # at the north pole c(0, ..., 0, 1)
  phi1 <- array(phi1, c(m, p, K))
  dimnames(phi1) <- list(t=pts, j=seq_len(p), k=seq_len(K))

  # Rotate to the locations of mu
  phi <- sapply(seq_along(pts), function(i) {
    R <- MakeRotMat(c(rep(0, p - 1), 1), mu[, i])
    R %*% phi1[i, , , drop=TRUE]
  }, simplify=FALSE)
  phi <- do.call(abind::abind, list(phi, along=0))
  dimnames(phi) <- dimnames(phi1)

  phi
}


MakePhi.Sphere <- 
  function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           type=c('cos', 'sin', 'fourier', 'legendre01')) {

  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }
  dimAmbient <- nrow(mu)
  dimIntrinsic <- calcIntDim.Sphere(mfd, dimAmbient)
  if (dimAmbient <= 1) {
    stop('mu must have more than 1 rows')
  }
  m <- length(pts)
  if (ncol(mu) != m) {
    stop('mu must have the same number of columns as the length of pts')
  }

  if (min(pts) < 0 || max(pts) > 1) {
    stop('The range of pts should be within [0, 1]')
  }

  ptsSqueeze <- do.call(c, lapply(seq_len(dimAmbient - 1), function(i) {
                          pts + (max(pts) + mean(diff(pts))) * (i - 1)
                        }))
  ptsSqueeze <- ptsSqueeze / max(ptsSqueeze)
  
  phi1 <- rbind(fdapace:::CreateBasis(K, ptsSqueeze, type),
                matrix(0, m, K)) / sqrt(dimIntrinsic) # at the north pole c(0, ..., 0, 1)
  phi1 <- array(phi1, c(m, dimAmbient, K))
  dimnames(phi1) <- list(t=pts, j=seq_len(dimAmbient), k=seq_len(K))

  # Rotate to the locations of mu
  basePt <- c(rep(0, dimIntrinsic), 1)
  phi <- sapply(seq_along(pts), function(i) {
    R <- MakeRotMat(basePt, mu[, i])
    R %*% phi1[i, , , drop=TRUE]
  }, simplify=FALSE)
  phi <- do.call(abind::abind, list(phi, along=0))
  dimnames(phi) <- dimnames(phi1)

  phi
}


MakePhi.SO <- 
  function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           type=c('cos', 'sin', 'fourier', 'legendre01')) {

  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }
  dimAmbient <- nrow(mu)
  dimTangent <- calcTanDim.SO(mfd, dimAmbient)
  dimMat <- round(sqrt(dimAmbient))
  if (dimAmbient <= 1) {
    stop('mu must have more than 1 rows')
  }
  m <- length(pts)
  if (ncol(mu) != m) {
    stop('mu must have the same number of columns as the length of pts')
  }

  if (min(pts) < 0 || max(pts) > 1) {
    stop('The range of pts should be within [0, 1]')
  }

  ptsSqueeze <- do.call(c, lapply(seq_len(dimTangent), function(i) {
                          pts + (max(pts) + mean(diff(pts))) * (i - 1)
                        }))
  ptsSqueeze <- ptsSqueeze / max(ptsSqueeze)
  
  phi <- array(fdapace:::CreateBasis(K, ptsSqueeze, type), c(m, dimTangent, K)) / 
    sqrt(dimTangent)

  dimnames(phi) <- list(t=pts, j=seq_len(dimTangent), k=seq_len(K))
  phi

  # vapply(seq_len(m), function(tt) {
    # vapply(seq_len(K), function(k) {
      # browser()
      # MakeSkewSym(epsFun)
    # }, rep(0, dimAmbient))
  # }, matrix(dimAmbient, K))
  # dimnames(phi1) <- list(t=pts, j=seq_len(dimAmbient), k=seq_len(K))

  # # Rotate to the locations of mu
  # basePt <- c(rep(0, dimTangent), 1)
  # phi <- sapply(seq_along(pts), function(i) {
    # R <- MakeRotMat(basePt, mu[, i])
    # R %*% phi1[i, , , drop=TRUE]
  # }, simplify=FALSE)
  # phi <- do.call(abind::abind, list(phi, along=0))
  # dimnames(phi) <- dimnames(phi1)

  # phi
}


## Generate a multivariate orthogonal basis on [0, 1].
# The indices of the returned array corresponds to (time, s, k)
MakePhi.Euclidean <- 
  function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           type=c('cos', 'sin', 'fourier', 'legendre01')) {

  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }
  dimAmbient <- nrow(mu)
  dimIntrinsic <- dimAmbient
  if (dimAmbient <= 1) {
    stop('mu must have more than 1 rows')
  }
  m <- length(pts)
  if (ncol(mu) != m) {
    stop('mu must have the same number of columns as the length of pts')
  }

  if (min(pts) < 0 || max(pts) > 1) {
    stop('The range of pts should be within [0, 1]')
  }

  ptsSqueeze <- do.call(c, lapply(seq_len(dimAmbient), function(i) {
                          pts + (max(pts) + mean(diff(pts))) * (i - 1)
                        }))
  ptsSqueeze <- ptsSqueeze / max(ptsSqueeze)
  
  phi1 <- fdapace:::CreateBasis(K, ptsSqueeze, type) / sqrt(dimAmbient)
  phi1 <- array(phi1, c(m, dimAmbient, K))
  dimnames(phi1) <- list(t=pts, j=seq_len(dimAmbient), k=seq_len(K))

  # Rotate to the locations of mu
  phi <- sapply(seq_along(pts), function(i) {
    phi1[i, , , drop=TRUE]
  }, simplify=FALSE)
  phi <- do.call(abind::abind, list(phi, along=0))
  dimnames(phi) <- dimnames(phi1)

  phi
}


## Generate an orthogonal basis on [0, 1] \times [0, 1], using the product basis functions \phi_j(t) \phi_l(s). Ordered by (j + l) and then j.
# The indices of the returned array corresponds to (time, s, k)
MakePhi.L2 <- 
  function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           type=c('cos', 'sin', 'fourier', 'legendre01')) {

  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }
  dimAmbient <- nrow(mu)
  dimTangent <- calcTanDim(mfd, dimAmbient)
  if (dimAmbient <= 1) {
    stop('mu must have more than 1 rows')
  }
  m <- length(pts)
  if (ncol(mu) != m) {
    stop('mu must have the same number of columns as the length of pts')
  }

  if (length(pts) != 1 && (min(pts) < 0 || max(pts) > 1)) {
    stop('The range of pts should be within [0, 1]')
  }

  ptss <- seq(0, 1, length.out=dimTangent)
  ptst <- pts

  if (m == 1) { # Scalar case
    phis <- fdapace:::CreateBasis(K, ptss, type)
    phi <- array(phis, c(m, dimTangent, K))
  } else {
    Keach <- ceiling(sqrt(2 * K))
    phis <- fdapace:::CreateBasis(Keach, ptss, type)
    phit <- fdapace:::CreateBasis(Keach, ptst, type)
    phi <- array(NA, c(m, dimTangent, K))
    k <- 0
    j <- 0
    l <- 1

    indm <- which(matrix(TRUE, Keach, Keach), arr.ind=TRUE)
    ord <- order(rowSums(indm))

    for (k in seq_len(K)) {
      ind <- indm[ord[k], ]
      j <- ind[1]
      l <- ind[2]
      phi[, , k] <- tcrossprod(phit[, j, drop=FALSE], phis[, l, drop=FALSE])
    }
  }
  dimnames(phi) <- list(t=pts, j=seq_len(dimAmbient), k=seq_len(K))

  phi
}


MakePhi.HS <- 
  function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           type=c('cos', 'sin', 'fourier', 'legendre01')) {

  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }
  dimAmbient <- nrow(mu)
  dimTangent <- calcTanDim(mfd, dimAmbient)
  if (dimAmbient <= 1) {
    stop('mu must have more than 1 rows')
  }
  m <- length(pts)
  if (ncol(mu) != m) {
    stop('mu must have the same number of columns as the length of pts')
  }

  if (length(pts) != 1 && (min(pts) < 0 || max(pts) > 1)) {
    stop('The range of pts should be within [0, 1]')
  }

  ptss <- seq(0, 1, length.out=dimTangent)
  ptst <- pts
  phi <- array(Inf, c(m, dimTangent, K))

  # Jt is the number of eigenfunctions corresponding to t
  # Jt <- 2
  # Js <- ceiling(K / Jt) + 1

  if (m == 1) { # Scalar case
    phis <- fdapace:::CreateBasis(K + 1, ptss, type)
    phi <- array(phis[, -1], c(m, dimTangent, K))
  } else {
    Jt <- Js <- ceiling(sqrt(2 * K))
    phis <- fdapace:::CreateBasis(Js, ptss, type)
    phit <- fdapace:::CreateBasis(Jt, ptst, type)

    # phi_1(t) phi1(s), phi_2(t) phi_1(s), phi_1(t) phi_2(s), etc
    matTRUE <- matrix(TRUE, Jt, Js) 
    indm <- which(matTRUE, arr.ind=TRUE)
    ord <- order(rowSums(indm))

    for (k in seq_len(K)) {
      # ind <- indm[k, ]
      ind <- indm[ord[k], ]
      j <- ind[1]
      l <- ind[2] + 1
      # browser()
      phi[, , k] <- tcrossprod(phit[, j, drop=FALSE], phis[, l, drop=FALSE])
    }
  }

  basePt <- Normalize(phis[, 1])

  dimnames(phi) <- list(t=pts, j=seq_len(dimAmbient), k=seq_len(K))

  # Rotate to the locations of mu
  res <- sapply(seq_along(pts), function(i) {
    R <- MakeRotMat(basePt, Normalize(mu[, i]))
    res <- R %*% phi[i, , , drop=TRUE]
    res1 <- projectTangent.HS(mfd, mu[, i], res)
    # if (mean(abs(res - res1)) > 1e-2) stop('Rot incorrect')
    # res2 <- projectTangent.HS(mfd, mu[, i], res[, 3])
    res1
  }, simplify=FALSE)
  res <- do.call(abind::abind, list(res, along=0))
  dimnames(res) <- dimnames(phi)

  res
}

# MakePhi.HS <- 
  # function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           # type=c('cos', 'sin', 'fourier', 'legendre01')) {

  # type <- match.arg(type)
  # if (!is.matrix(mu)) {
    # stop('mu has to be a matrix')
  # }
  # dimAmbient <- nrow(mu)
  # dimTangent <- calcTanDim(mfd, dimAmbient)
  # if (dimAmbient <= 1) {
    # stop('mu must have more than 1 rows')
  # }
  # m <- length(pts)
  # if (ncol(mu) != m) {
    # stop('mu must have the same number of columns as the length of pts')
  # }

  # if (min(pts) < 0 || max(pts) > 1) {
    # stop('The range of pts should be within [0, 1]')
  # }

  # ptss <- seq(0, 1, length.out=dimTangent)
  # ptst <- pts
  # Keach <- ceiling(sqrt(2 * K))
  # phis <- fdapace:::CreateBasis(Keach, ptss, type)
  # phit <- fdapace:::CreateBasis(Keach, ptst, type)
  # phi <- array(NA, c(m, dimTangent, K))

  # matTRUE <- matrix(TRUE, Keach, Keach)

  # indm <- which(matTRUE, arr.ind=TRUE)
  # ord <- order(rowSums(indm))

  # # Ensure the bases are orthogonal to basePt. Take only phi_2, ..., phi_{K+1},
  # # and then rotate them from phi_1 to mu
  # if (type %in% c('fourier', 'cos', 'legendre01', 'sin')) {
    # lplus <- 1
  # } # else {
    # # lplus <- 0
  # # }

  # for (k in seq_len(K)) {
    # ind <- indm[ord[k], ]
    # j <- ind[1]
    # l <- ind[2] + lplus
    # phi[, , k] <- tcrossprod(phit[, j, drop=FALSE], phis[, l, drop=FALSE])
    # # if (k == 3) browser()
  # }

  # basePt <- Normalize(phis[, 1])

  # dimnames(phi) <- list(t=pts, j=seq_len(dimAmbient), k=seq_len(K))

  # # Rotate to the locations of mu
  # res <- sapply(seq_along(pts), function(i) {
    # R <- MakeRotMat(basePt, Normalize(mu[, i]))
    # res <- R %*% phi[i, , , drop=TRUE]
    # res1 <- projectTangent.HS(mfd, mu[, i], res)
    # # if (max(abs(res - res1)) > 1e-2) browser()
    # # res2 <- projectTangent.HS(mfd, mu[, i], res[, 3])
    # res1
  # }, simplify=FALSE)
  # res <- do.call(abind::abind, list(res, along=0))
  # dimnames(res) <- dimnames(phi)

  # res
# }


# The dimensions of phi corresponds to (time, j, k)
MakePhi <- function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           type=c('cos', 'sin', 'fourier', 'legendre01')) {
  UseMethod('MakePhi', mfd)
}


# For each time point (column in mu), generate n isotropic noises distributed as epsFun on the tangent space at mu(t)
# return an array with size (n, dimIntrinsic, nPts) 
NoiseTangent <- function(mfd, n, mu, sigma2=0, epsFun = rnorm) {
  UseMethod('NoiseTangent', mfd)
}


NoiseTangent.Euclidean <- function(mfd, n, mu, sigma2=0, epsFun = rnorm) {

  m <- ncol(mu)
  dimIntrinsic <- dimAmbient <- nrow(mu)
  res <- vapply(seq_len(m), function(tt) {
           eps <- matrix(epsFun(n * dimIntrinsic), n, dimIntrinsic)
           eps
         }, matrix(0, n, dimAmbient)) * sqrt(sigma2)
  # res <- aperm(res, c(1, 3, 2))
  res
}


NoiseTangent.Sphere <- function(mfd, n, mu, sigma2=0, epsFun = rnorm) {

  m <- ncol(mu)
  dimAmbient <- nrow(mu)
  res <- vapply(seq_len(m), function(tt) {
           rot <- MakeRotMat(c(rep(0, dimAmbient - 1), 1), mu[, tt])
           eps <- cbind(matrix(epsFun(n * (dimAmbient - 1)), n, dimAmbient - 1), 0) %*% t(rot)
           eps
         }, matrix(0, n, dimAmbient)) * sqrt(sigma2)
  # res <- aperm(res, c(1, 3, 2))
  res
}


NoiseTangent.SO <- function(mfd, n, mu, sigma2=0, epsFun = rnorm) {

  m <- ncol(mu)
  dimAmbient <- nrow(mu)
  dimTangent <- calcTanDim.SO(mfd, dimAmbient)
  res <- array(epsFun(n * dimTangent * m), c(n, dimTangent, m)) * sqrt(sigma2)
  res

}


NoiseTangent.SO <- function(mfd, n, mu, sigma2=0, epsFun = rnorm) {

  m <- ncol(mu)
  dimAmbient <- nrow(mu)
  dimTangent <- calcTanDim.SO(mfd, dimAmbient)
  res <- array(epsFun(n * dimTangent * m), c(n, dimTangent, m)) * sqrt(sigma2)
  res

}





#' Simulate a manifold-valued process
#'
#' @param mfd An object whose class specifies the manifold to use
#' @param n The number of curves
#' @param mu A matrix that specifies the mean function. Each column stands for a manifold-valued data point
#' @param pts The time points
#' @param K The number of components to simulate
#' @param lambda The eigenvalues
#' @param sigma2 The variance for the isotropic measurement errors
#' @param xiFun The distribution of the RFPC scores
#' @param epsFun The distribution of the measurement errors
#' @param epsBase The base point of the (tangent) measurement errors.
#' @param basisType The type of basis on the tangent spaces
#' @export
MakeMfdProcess <- function(mfd, 
  n, mu, pts=seq(0, 1, length.out=ncol(mu)), K=2, lambda=rep(1, K), sigma2=0, 
  xiFun = rnorm, epsFun = rnorm, epsBase = c('mu', 'X'), 
  basisType=c('cos', 'sin', 'fourier', 'legendre01')) {

  epsBase <- match.arg(epsBase)
  basisType <- match.arg(basisType)
  dimAmbient <- nrow(mu)
  dimTangent <- calcTanDim(mfd, dimAmbient)
  m <- length(pts)

  # A function that generates iid centered and scaled random variables for scores
  xi <- matrix(xiFun(n * K), n, K) %*% diag(sqrt(lambda), nrow = K)
  phi <- matrix(MakePhi(mfd, K, mu, pts, basisType), ncol=K)
  xiphi <- aperm(array(xi %*% t(phi), c(n, m, dimTangent)), c(1, 3, 2))

  X <- sapply(seq_len(m), function(tt) {
    mu0 <- mu[, tt]
    V <- matrix(xiphi[, , tt], n, dimTangent)
    t(rieExp(mfd, mu0, t(V)))
  }, simplify='array')
  dimnames(X) <- list(i = seq_len(n), j = seq_len(dimAmbient), t = pts)

  if (sigma2 > 1e-15) {
    if (epsBase == 'mu') {
      ## Rewrite here
      epsArray <- NoiseTangent(mfd, n, mu, sigma2, epsFun)
      xiphiN <- xiphi + epsArray
      XNoisy <- sapply(seq_len(m), function(tt) {
        mu0 <- mu[, tt]
        V <- t(matrix(xiphiN[, , tt], n, dimTangent))
        t(rieExp(mfd, mu0, V))
      }, simplify='array')
    } else if (epsBase == 'X') { # Does not seem to be useful
      stop('not implemented now')
      # XNoisy <- apply(X, c(1, 3), function(x) {
        # rot <- MakeRotMat(c(rep(0, dimTangent - 1), 1), x)
        # eps <- rot %*% matrix(c(epsFun(dimTangent - 1), 0) * sqrt(sigma2))
        # t(rieExp(mfd, x, eps))
      # })
      # XNoisy <- aperm(XNoisy, c(2, 1, 3))
    }
  } else {
    XNoisy <- X
  }

  res <- list(XNoisy = XNoisy, X = X, T = pts, xi = xi, phi=phi)
  res
}



## Generate process on a sphere by exponential maps
# Return a three dimensional array X[i, j, t] = X_{ij}(t)
# sigma2: The entriwise (isotropic) variance for added errors
# xiFun: A function to generate xi. The first argument is the number of samples to generate.
# epsFun: A function to generate noise, similar to xiFun.
# epsBase: What is the base point p for adding noises. The noisy observations are Exp_p(eps), where p = 'mu' or 'X'
MakeSphericalProcess <- function(
  n, mu, pts=seq(0, 1, length.out=ncol(mu)), K=2, lambda=rep(1, K), sigma2=0, 
  xiFun = rnorm, epsFun = rnorm, epsBase = c('mu', 'X'), 
  basisType=c('cos', 'sin', 'fourier', 'legendre01')) {

  epsBase <- match.arg(epsBase)
  basisType <- match.arg(basisType)
  p <- nrow(mu)
  m <- length(pts)

  # A function that generates iid centered and scaled random variables for scores
  xi <- matrix(xiFun(n * K), n, K) %*% diag(sqrt(lambda), nrow = K)
  phi <- matrix(MakephiS(K, mu, pts, basisType), ncol=K)
  xiphi <- array(xi %*% t(phi), c(n, m, p))

  mfd <- structure(1, class='Sphere')
  X <- sapply(seq_len(m), function(tt) {
    mu0 <- mu[, tt]
    V <- matrix(xiphi[, tt, ], dim(xiphi)[1], dim(xiphi)[3])
    t(rieExp(mfd, mu0, t(V)))
  }, simplify='array')
  dimnames(X) <- list(i = seq_len(n), j = seq_len(p), t = pts)

  if (sigma2 > 1e-15) {
    if (epsBase == 'mu') {
      epsArray <- vapply(seq_len(m), function(tt) {
        rot <- MakeRotMat(c(rep(0, p - 1), 1), mu[, tt])
        eps <- cbind(matrix(epsFun(n * (p - 1)), n, p - 1), 0) %*% t(rot) * 
          sqrt(sigma2)
        eps
      }, matrix(0, n, p))
      xiphiN <- xiphi + aperm(epsArray, c(1, 3, 2))
      XNoisy <- sapply(seq_len(m), function(tt) {
        mu0 <- mu[, tt]
        V <- matrix(xiphiN[, tt, ], dim(xiphiN)[1], dim(xiphiN)[3])
        t(rieExp(mfd, mu0, t(V)))
      }, simplify='array')
    } else if (epsBase == 'X') {
      XNoisy <- apply(X, c(1, 3), function(x) {
        rot <- MakeRotMat(c(rep(0, p - 1), 1), x)
        eps <- rot %*% matrix(c(epsFun(p - 1), 0) * sqrt(sigma2))
        t(rieExp(mfd, x, eps))
      })
      XNoisy <- aperm(XNoisy, c(2, 1, 3))
    }
  } else {
    XNoisy <- X
  }


  res <- list(XNoisy = XNoisy, X = X, T = pts, xi = xi, phi=phi)
  res
}


#' Sparsify a dense Riemannian process
#' @param X A three dimensional array such that X[i, j, t] = X_{ij}(t). 
#' @param pts A vector of grid points corresponding to the columns of \code{samp}.
#' @param sparsity A vector of integers. The number of observation per sample is chosen to be one of the elements in sparsity with equal chance
#' @export
SparsifyM <- function(X, pts, sparsity) {

  if (!is.array(X) || length(dim(X)) != 3) {
      stop("X needs to be a 3D array")
  }

  if (dim(X)[3] != length(pts)) {
      stop("The number of time points in X needs to be equal to the length of pts")
  }

  if (length(sparsity) == 1) {
      sparsity <- c(sparsity, sparsity)
  }

  indEach <- lapply(seq_len(dim(X)[1]), function(x) sort(sample(dim(X)[3], 
    sample(sparsity, 1))))

  Lt <- lapply(indEach, function(x) pts[x])
  Ly <- lapply(seq_along(indEach), function(x) {
      ind <- indEach[[x]]
      y <- matrix(X[x, , ind], ncol=length(ind))
      unname(y)
  })

  return(list(Lt = Lt, Ly = Ly))
}


CFPCA <- function(Ly, Lt, optns=list()) {

  # Componentwise FPCA and then summarize
  p <- nrow(Ly[[1]])

  KUse <- optns[['KUse']]
  optns <- optns[names(optns) != 'KUse']
  if (is.null(KUse)) {
    KUse <- max(sapply(resFPCA, function(res) length(res$lambda)))
  }

  resFPCA <- lapply(seq_len(p), function(j) {
    Ly <- lapply(Ly, `[`, j, )
    res <- FPCA(Ly, Lt, optns)
    res
  })

  mfd <- structure(1, class='Euclidean')
  workGrid <- resFPCA[[1]]$workGrid
  obsGrid <- resFPCA[[1]]$obsGrid
  muEst <- t(vapply(resFPCA, function(res) {
    muObs <- ConvertSupport(workGrid, obsGrid, res[['mu']])
  }, obsGrid))

  xiFPCA <- do.call(cbind, lapply(resFPCA, `[[`, 'xiEst'))
  if (ncol(xiFPCA) < KUse) {
    stop('Cannot obtain KUse components')
  }
  xiPCA <- prcomp(xiFPCA, center=FALSE)
  xiEst <- xiPCA[['x']][, seq_len(KUse), drop=FALSE] #%*% xiPCA[['rotation']][, seq_len(KUse), drop=FALSE]
  lam <- xiPCA[['sdev']][seq_len(KUse)]

  phiFPCA <- do.call(cbind, lapply(resFPCA, function(res) {
    apply(res$phi, 2, function(x) {
            ConvertSupport(workGrid, obsGrid, mu=x)
    })
  }))
  KEach <- sapply(resFPCA, function(res) length(res$lambda))
  indEach <- split(seq_len(sum(KEach)), rep(seq_len(p), KEach))

  phiEst <- vapply(indEach, function(ind) {
    phi <- phiFPCA[, ind, drop=FALSE]
    rot <- xiPCA[['rotation']][ind, seq_len(KUse), drop=FALSE]
    phi %*% rot
  }, matrix(0, nrow(phiFPCA), KUse))
  phiEst <- aperm(phiEst, c(1, 3, 2))

  res <- list(muObs = muEst, 
              phiObs = phiEst, 
              lam = lam, 
              xi=xiEst,
              K = KUse, 
              mfd = mfd, 
              optns=optns)

  res
}


# Calculate errors
# phiTrue, phiEst: matrices, whose columns are eigenfunctions
ErrPhiXi <- function(phiTrue, xiTrue, phiEst, xiEst) {

  if (!is.matrix(phiTrue)) phiTrue <- matrix(phiTrue)
  if (!is.matrix(xiTrue)) xiTrue <- matrix(xiTrue)
  if (!is.matrix(phiEst)) phiEst <- matrix(phiEst)
  if (!is.matrix(xiEst)) xiEst <- matrix(xiEst)

  flip <- diag(crossprod(phiTrue, phiEst)) < 0
  phiEst[, flip] <- -phiEst[, flip]
  xiEst[, flip] <- -xiEst[, flip]

  errPhi <- colMeans((phiTrue - phiEst)^2)
  errXi <- colMeans((xiTrue - xiEst)^2)

  list(phi=errPhi, xi=errXi, flip=flip)
}


# xTrue: A 3-D array, whose dimensions corresponds to (i, j, t)
# phiEst: A 3-D array, whose dimensions corresponds to (t, j, k), where k is the index for eigenfunction.
# RFPCA: If true, fit values on the manifold through the exponential maps; if false, fit in the ambient space
# dist: 'distance', just use the manifold distance; 'sqDist', square the data and fitted values, and then apply L2 distance
# mfd The manifold on which the geodesic distance should be calculated or should be projected to 
ErrFitted <- function(mfd, xTrue, object, newLy, newLt, RFPCA=TRUE, 
                      expMap=TRUE, dist=c('distance', 'sqDist'), Kvec, muTrue, rel=FALSE) {

  dist <- match.arg(dist)
  if (!missing(RFPCA)) {
    expMap <- RFPCA
    warning('Option `RFPCA` is deprecated. Use `expMap` instead.')
  }

  if (!expMap) {
    if (inherits(object[['mfd']], 'L2')) {
      object[['mfd']] <- structure(1, class='L2')
    } else {
      object[['mfd']] <- structure(1, class='Euclidean')
    }
  }

  if (!missing(newLt) && !missing(newLy)) { # predict
    fitFunc <- function(K) predict.RFPCA(object, newLy, newLt, K=K, xiMethod='IN', type='traj')
  } else { # fit
    fitFunc <- function(K) fitted.RFPCA(object, K=K, grid='obs')
  }

  # if (!is.matrix(phiEst)) phiEst <- matrix(phiEst)
  # if (!is.matrix(xiEst)) xiEst <- matrix(xiEst)

  # if (length(dim(xTrue)) > 2 || length(dim(phiEst)) > 2) {
    # stop('Array valued input are not supported')
  # }

  n <- dim(xTrue)[1]
  p <- dim(xTrue)[2]
  m <- dim(xTrue)[3]

  if (missing(Kvec)) {
    Kvec <- seq_len(ncol(object[['xi']]))
  }

  if (dist == 'distance') {
    dfunc <- function(a, b) distance(mfd, a, b)
  } else if (dist == 'sqDist') {
    mfdL2 <- structure(1, class='L2')
    dfunc <- function(a, b) distance(mfdL2, a^2, b^2)
  }

  res <- sapply(Kvec, function(K) {

    fit <- fitFunc(K)
    if (!expMap) {
      fit <- aperm(apply(fit, c(1, 3), project, mfd=mfd), c(2, 1, 3))
    }

      # if (pred) { # Get predicted values
        # if (!expMap) {
          # xiphi <- tcrossprod(xiEst[, seq_len(K), drop=FALSE], 
                              # matrix(phiEst[, , seq_len(K)], m * p, K))
          # fit <- matrix(t(muEst), n, p * m, byrow=TRUE) + xiphi
          # fit <- aperm(apply(array(fit, c(n, m, p)), c(1, 2), project, mfd=mfd), c(2, 1, 3))
        # } else {
          # fit <- predict(mfd, muEst, xiEst, phiEst, K)
        # }
      # } else { # Get fitted values
        # if (!expMap) {
          # xiphi <- tcrossprod(xiEst[, seq_len(K), drop=FALSE], 
                              # matrix(phiEst[, , seq_len(K)], m * p, K))
          # fit <- matrix(t(muEst), n, p * m, byrow=TRUE) + xiphi
          # fit <- aperm(apply(array(fit, c(n, m, p)), c(1, 2), project, mfd=mfd), c(2, 1, 3))
        # } else {
          # # mfd <- structure(1, class='Sphere')
          # fit <- Fitted(mfd, muEst, xiEst, phiEst, K)
        # }
      # }
    resi <- sapply(seq_len(n), function(i) {
                  mean(dfunc(matrix(xTrue[i, , ], p, m), 
                             matrix(fit[i, , ], p, m))^2)
               })
    if (rel) {
      distFromMu <- sapply(seq_len(n), function(i) {
                             mean(dfunc(matrix(xTrue[i, , ], p, m), 
                                        muTrue)^2)
                           })
      resi <- resi / distFromMu
    }
    mean(resi)
  })

  res
}


# # pts: the support of xTrue
# ErrFittedFPCA <- function(mfd, xTrue, pts, resFPCA, KUse) {

  # n <- dim(xTrue)[1]
  # p <- dim(xTrue)[2]
  # m <- dim(xTrue)[3]

  # if (missing(KUse)) {
    # KUse <- seq_len(max(sapply(resFPCA, function(res) length(res$lambda))))
  # }

  # muEst <- t(apply(sapply(resFPCA, `[[`, 'mu'), 2, function(x) {
                  # ConvertSupport(resFPCA[[1]]$workGrid, pts, mu=x)
  # }))

  # xiFPCA <- do.call(cbind, lapply(resFPCA, `[[`, 'xiEst'))
  # xiPCA <- prcomp(xiFPCA)
  # xiEst <- xiPCA$x[, seq_len(KUse), drop=FALSE]
  # xiEst <- sapply(KUse, function(K) {
    # tcrossprod(xiPCA$x[, seq_len(K), drop=FALSE], xiPCA$rotation[, seq_len(K), drop=FALSE])
  # }, simplify='array')

  # phiFPCA <- lapply(resFPCA, function(res) {
    # apply(res$phi, 2, function(x) {
            # ConvertSupport(res$workGrid, pts, mu=x)
    # })
  # })

  # KEach <- sapply(resFPCA, function(res) length(res$lambda))
  # indEach <- split(seq_len(sum(KEach)), rep(seq_len(p), KEach))

  # res <- sapply(KUse, function(K) {

    # xi <- lapply(indEach, function(ind) xiEst[, ind, K])
    # xiphi <- sapply(seq_len(p), function(j) {
                      # tcrossprod(xi[[j]], phiFPCA[[j]])
    # }, simplify='array')
    # xiphi <- aperm(xiphi, c(1, 3, 2))

      # fit <- array(apply(xiphi, 1, function(mati) {
                     # res <- muEst + mati
                     # # if (mapToSphere) {
                       # res <- apply(res, 2, Normalize)
                     # # }
                     # res
                   # }), c(p, m, n))
      # fit <- aperm(fit, c(3, 1, 2))

    # mean(sapply(seq_len(n), function(i) {
                  # mean(distance(mfd, 
                                # matrix(xTrue[i, , ], p, m), 
                                # matrix(fit[i, , ], p, m))^2
                 # )
               # }))
    # })

  # res
# }


ErrMu <- function(mfd, muTrue, muEst) {
  muEst <- apply(muEst, 2, project, mfd=mfd)
  mean(distance(mfd, muTrue, muEst)^2)
}


GetRVGenerator <- function(name='norm') {

  if (name == 'norm') {
    rnorm
  } else if (name == 'exp') {
    function(n) rexp(n) - 1
  } else if (name == 't10') {
    function(n) rt(n, df=10) * sqrt(8 / 10) 
  } else {
    stop(sprintf("'%s' random number generator not implemented", name))
  }

}


GetSettingName <- function(settings, digits=3, display=c('short', 'pretty')) {

  display <- match.arg(display)

  if (is.null(names(settings))) {
    stop("'settings' must be a named list")
  }

  type <- sapply(settings, function(x) {
    if (is.list(x) && length(x) == 1) {
      x <- unlist(x)
    }

    switch(class(x),
      'character' = '%s', 
      'numeric' = sprintf('%%.%dg', digits), 
      'integer' = '%d', 
      stop('Not supported')
      )
  })
  val <- sapply(settings, function(x) {
    if (is.list(x) && length(x) == 1) {
      x <- unlist(x)
    }

    switch(class(x),
      'function' = stop('Not supported'), 
      'factor' = stop('Not supported'), 
      'integer' = round(mean(x)), 
      x[1]
      )
  }, simplify=FALSE)

  if (display == 'short') {
    form <- paste(paste0(names(settings), type), collapse='_')
  } else if (display == 'pretty') {
    form <- paste(paste0(names(settings), '=', type), collapse=', ')
  }

  fn <- do.call(sprintf, c(list(form), val))
  fn
}


## Old functions
## Generate eigenfunctions V(t) from a certain basis around a mean curve. 
# mu needs to correspond to pts.
MakephiSO <- function(K, mu, pts=seq(0, 1, length.out=nrow(mu)), 
                      type=c('cos', 'sin', 'fourier', 'legendre01')) {

  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }

  p <- ncol(mu) # ambient dimension
  d <- round(sqrt(p))
  if (p <= 1) {
    stop('mu must have more than 1 columns')
  }
  stopifnot(d == sqrt(p))

  p0 <- d * (d - 1) / 2 # true dimension

  m <- length(pts)
  if (nrow(mu) != m) {
    stop('mu must have the same number of rows as the length of pts')
  }

  if (min(pts) < 0 || max(pts) > 1) {
    stop('The range of pts should be within [0, 1]')
  }

  ptsSqueeze <- do.call(c, lapply(seq_len(p0), function(i) {
                          pts + (max(pts) + mean(diff(pts))) * (i - 1)
                        }))
  ptsSqueeze <- ptsSqueeze / max(ptsSqueeze)
  
  phip0 <- array(fdapace:::CreateBasis(K, ptsSqueeze, type), c(m, p0, K)) / sqrt(p0)

  phi <- sapply(seq_along(pts), function(i) {
                  mat <- matrix(phip0[i, , ], p0, K)
                  # muMat <- matrix(mu[i, ], d, d)

                  apply(mat, 2, function(v) {
                    MakeSkewSym(v)
                  })
                }, simplify=FALSE)
  phi <- do.call(abind::abind, list(phi, along=0)) / sqrt(2) # due to symmetrization
  dimnames(phi) <- list(t=pts, j=seq_len(p), k=seq_len(K))

  phi
}


# ## Generate mean function by exponential map. 
# # VtList: a list of coordinate functions for V(t). For type = 'SO' only the lower triangle should be specified.
# # p0: the base point of the exponential map
# Makemu.SO <- function(mfd, VtList, p0, pts=seq(0, 1, length.out=50)) {

  # Vmat <- sapply(VtList, function(f) f(pts)) # lower triangle only
  # res <- rieExp(mfd, p0, t(Vmat))

  # res
# }


# Project an RFPCA object to Hilbert Sphere
# Does not work if truncated
ProjectObjectHS <- function(obj) {
  mfdHS <- structure(1, class='HS')
  obj$mfdHS <- mfdHS
  obj$muReg <- project(mfdHS, obj$muReg)
  obj$muWork <- project(mfdHS, obj$muWork)
  obj$muObs <- project(mfdHS, obj$muObs)
  obj$cov <- projectCov(mfdHS, obj$cov, obj$muWork)
  obj$covObs <- projectCov(mfdHS, obj$covObs, obj$muObs)
  obj$phiObsTrunc <- projectPhi(mfdHS, obj$phiObsTrunc, obj$muObs)

  obj
}


HSToRootDens <- function(x) {

  p <- length(x)

  if (all(x < 0)) {
    warning('A negative function cannot be projected to a density. Return the uniform density instead')
    x <- rep(sqrt(p - 1) / sqrt(p), p)
  } else {
    x[x < 0] <- 0
    x <- x * sqrt((p - 1) / sum(x^2))
  }

  x
}


HSSampToDensities <- function(samp) {

  # n <- dim(samp)[1]
  # p <- dim(samp)[2]
  # m <- dim(samp)[3]
  dn <- dimnames(samp$X)
  samp$X <- aperm(apply(samp$X, c(1, 3), HSToRootDens), c(2, 1, 3))
  samp$XNoisy <- aperm(apply(samp$XNoisy, c(1, 3), HSToRootDens), c(2, 1, 3))
  dimnames(samp$X) <- dimnames(samp$XNoisy) <- dn
  samp

}
