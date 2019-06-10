## Class for general SO(n) manifolds. In this implementation, we follow the convention that
## we use only the lower triangle of the skew-symmetric matrices on the tangent spaces to 
## denote the tangent vectors. 

# #' Riemannian metric of tangent vectors
# #' 
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param p The base point which could be NULL if the metric does not rely on the base points
# #' @param U A D*n matrix, each column represents a tangent vector at \emph{p}
# #' @param V A D*n matrix, each column represents a tangent vector at \emph{p}
# #' 
# #' @details The tangent vectors can be represented in a coordinate frame, or in ambient space
# #' @export
metric.SO <- function(mfd,U,V,p=NULL) {

  U <- as.matrix(U)
  V <- as.matrix(V)
  res <- colSums(U * V)
  return(res)
}

# #' The norm induced by the Riemannian metric tensor on tangent spaces
# #' 
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param p The base point which could be NULL if the norm does not rely on it
# #' @param U A D*n matrix, where each column represents a tangent vector at \emph{p}
norm.SO <- function(mfd,U,p=NULL) {
    return(sqrt(colSums(U*U)))
}


# #' Geodesic distance of points on the manifold
# #' 
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param X A D*n matrix. Each column represents a point on manifold
# #' @param Y A D*n matrix. Each column represents a point on manifold
# #' @return A 1*n vector. The \emph{i}th element is the geodesic distance of \code{X[,i]} and \code{Y[,i]}
distance.SO <- function(mfd, X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  d <- round(sqrt(nrow(X)))
  n <- ncol(X)

  V <- sapply(seq_len(n), function(i) {
    x <- matrix(X[, i], d, d)
    y <- matrix(Y[, i], d, d)
    log.SO(X=as.numeric(crossprod(x, y)))
  })

  norm.SO(mfd, V)
}


# #' Geodesic curve stating at a point
# #' 
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param p The starting point of the geodesic curve
# #' @param h A matrix with each column representing a tangent vector. If there is only one tangent vector is supplied, then it is replicated to match the length of \emph{t}
# #' @param t A vector, the time points where the curve is evaluated.
# #' 
# #' @details The curve is \eqn{\gamma(t)=\mathrm{Exp}_p(th)}
# #' 
# #' @return A matrix with each colum representing a point on the manifold

# geodesicCurve.SO <- function(mfd,p,h,t) {
    # if(!is.matrix(p)) p <- as.matrix(p)
    # if(!is.matrix(h)) h <- as.matrix(h)
    # stopifnot(all(abs(crossprod(p, h) - 0) < 1e-8))
    # n <- dim(h)[2]
    # d <- dim(h)[1]
    # m <- length(t)

    # if (n==1) {
      # h <- h %*% matrix(t, nrow=1)
    # } else if (m != n) {
      # stop('If there is more than one tangent vectors, the number of tangent vectors must be equal to the number of time points')
    # }

    # exp.SO(mfd, p, h)
    
# }


# #' Riemannian exponential map at a point
exp.SO <- function(mfd, p, V, tol=1e-10) {

  if (!is.matrix(V)) {
    V <- matrix(V, ncol=1)
  }

  n <- ncol(V)
  d <- calcGeomPar.SO(dimTangent=nrow(V))

  if (missing(p)) {
    p <- matrix(diag(d), d^2, ncol(V))
  } else {
    if (is.matrix(p)) {
      stopifnot(nrow(p) == d^2)
    } else {
      p <- matrix(p)
    }
    if (ncol(p) == 1) {
      p <- matrix(p, nrow(p), n)
    }
  }

  # each col of V needs to be orthogonal to p

  res <- sapply(seq_len(n), function(i) {
                  vMat <- MakeSkewSym(V[, i])
                  muMat <- matrix(p[, i], d, d)
                  stopifnot(isSkewSym(vMat))

                  if (d == 3) {
                    ExpMSO3(vMat, tol) %*% muMat
                  } else {
                    # seems only 'Eigen' works but not 'Higham08'
                    ExpM(vMat, method='R_Eigen') %*% muMat
                  }
               })

  res
}


# #' Riemannian log map at a point
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param p A matrix. Each column represents a base point of the log map. If only one base point is supplied, then it is replicated to match the number of points in \emph{X}. If missing, the basepoint is set to be the identity matrix.
# #' @param X A matrix. Each column represents a point on the manifold
# #' @return A matrix with the \emph{i}th column being the log map of the \emph{i}th point

log.SO <- function(mfd, p, X, tol=1e-10) {

  if (!is.matrix(X)) {
    X <- matrix(X, ncol=1)
  }

  n <- ncol(X)
  d <- round(sqrt(nrow(X)))

  if (missing(p)) {
    p <- matrix(diag(d), nrow(X), ncol(X))
  } else {
    if (is.matrix(p)) {
      stopifnot(nrow(p) == nrow(X))
    } else {
      p <- matrix(p)
    }
    if (ncol(p) == 1) {
      p <- matrix(p, nrow(X), n)
    }
  }

  res <- sapply(seq_len(n), function(i) {
                  xMat <- matrix(X[, i], d, d)
                  muMat <- matrix(p[, i], d, d)
                  stopifnot(isSO(xMat, tol))
                  if (d == 3) {
                    LogMSO3(tcrossprod(xMat, muMat))
                  } else {
                    # seems only 'Eigen' works but not 'Higham08'
                    LogM(tcrossprod(xMat, muMat), method='Eigen')
                  }
                 })

  res[c(lower.tri(diag(d))), , drop=FALSE]
}

# Matrix logarithm
# mat A matrix
LogM <- function(mat, method='Higham08') {

  expm::logm(mat, method)

}

# Works only for d = 3
LogMSO3 <- function(mat, tol=1e-10) {

  if (nrow(mat) == 3) {
    theta <- acos((sum(diag(mat)) - 1) / 2)
    if (abs(theta) < tol) {
      res <- matrix(0, 3, 3)
    } else {
      res <- theta / (2 * sin(theta)) * (mat - t(mat))
    }
  } else {
    stop()
  }
  res
}


ExpM <- function(mat, method='R_Eigen') {
  expm::expm(mat, method)
}


ExpMSO3 <- function(mat, tol=1e-10) {

  if (nrow(mat) == 3) {
    a <- sqrt(0.5 * sum(mat^2))
    if (abs(a) < tol) {
      res <- diag(3)
    } else {
      res <- diag(3) + sin(a) / a * mat + (1 - cos(a)) / a^2 * mat %*% mat
    }
  } else {
    stop()
  }
  res
}


# Project extrinsic tangent space data onto the tangent space at p
# p must be a single point on the manifold
projectTangent.SO <- function(mfd, p, X, projMatOnly=FALSE) {

  d <- round(sqrt(length(p)))
  stopifnot(isSO(matrix(p, d, d)))
  projMat <- diag(d)

  if (projMatOnly) {
    return(projMat)
  } else {
    return(X)
  }
}


## Helper functions from sphere/R/func.R
## Functions for rotational groups SO(n). c.f. Rahman et al 2006

isOrth <- function(x, tol=1e-14) {
  if (length(x) == 1) {
    x <- matrix(x, 1, 1)
  }

  d <- nrow(x)
  max(abs(crossprod(x) - diag(d))) <= tol
}


isSO <- function(x, tol=1e-14) {
  if (length(x) == 1) {
    x <- matrix(x, 1, 1)
  }

  d <- nrow(x)
  max(abs(crossprod(x) - diag(d))) <= tol && abs(det(x) - 1) <= tol
}


isSkewSym <- function(x, tol=1e-14) {
  max(abs(x + t(x) - 0)) <= tol
}


# Make a skew-symmetric matrix by specifying the lower trianglular elements.
MakeSkewSym <- function(v) {

  p <- length(v)
  d <- calcGeomPar.SO(dimTangent=p)

  A <- matrix(0, d, d)
  A[lower.tri(A)] <- v
  A <- A - t(A)

  A
}


# Make a skew-symmetric matrix of functions by specifying the lower trianglular functions.
MakeSkewSymFunc <- function(fList) {

    p <- length(fList)
    d <- calcGeomPar.SO(dimTangent=p)

    outList <- matrix(list(), d, d)
    outList[] <- lapply(outList, function(l) {function(x) rep(0, length(x))})
    outList[lower.tri(outList)] <- lapply(fList, function(f) {
                                          function(x) -1 * f(x)
                                        }) 
    outList <- t(outList)
    outList[lower.tri(outList)] <- fList

    outList

}


# Make an orthogonal matrix by specifying the lower trianglular elements of the skew-symmetric logarithm matrix.  
MakeSOMat <- function(v, method='Higham08.b') {

  expm::expm(MakeSkewSym(v), method)

}


DistSO <- function(x1, x2) {

  if (!is.matrix(x1)) {
    d1 <- round(sqrt(length(x1)))
    x1 <- matrix(x1, d1, d1)
  }
  if (!is.matrix(x2)) {
    d2 <- round(sqrt(length(x2)))
    x2 <- matrix(x2, d2, d2)
  }

  if (nrow(x1) == 3) {
    R <- crossprod(x1, x2)
    x <- (sum(diag(R)) - 1) / 2
    if (x > 1) {
      x <- 1
    } else if (x < -1) {
      x <- -1
    }
    sqrt(2) * abs(acos(x))
  } else {
    drop(sqrt(crossprod(as.numeric(expm::logm(crossprod(x1, x2), method='Eigen')))))
  }

}


# ExpMapSO <- function(V, mu, ...) {

  # if (is.matrix(mu)) {
    # stopifnot(dim(mu)[1] == dim(mu)[2])
    # d <- nrow(mu)
  # } else {
    # d <- round(sqrt(length(mu)))
    # stopifnot(identical(d, sqrt(length(mu))))
  # }

  # if (!is.matrix(V)) {
    # V <- matrix(V, nrow=1)
  # }

  # stopifnot(ncol(V) == length(mu))

  # muMat <- matrix(mu, d, d)
  # stopifnot(all.equal(crossprod(muMat), diag(d)))

  # res <- t(apply(V, 1, function(v) {
                   # vMat <- matrix(v, d, d)
                   # expm::expm(vMat, ...) %*% muMat
                 # }))

  # res
# }


# ExpMap1 <- function(v, mu, type, ...) {

  # # stopifnot(length(v) == length(mu))

  # if (type == 'SO') {
    # d <- round(sqrt(length(mu)))

    # vMat <- matrix(v, d, d)
    # muMat <- matrix(mu, d, d)
    # res <- as.numeric(expm::expm(vMat, ...) %*% muMat)
  # }

  # res
# }


# # works only for SO(3). c.f. Rahman et al 2006 and https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
# ExpMap2 <- function(v, mu, type, ...) {

  # # stopifnot(length(v) == length(mu))
  # # stopifnot(length(v) == 9)

  # if (type == 'SO') {
    # d <- round(sqrt(length(mu)))

    # vMat <- matrix(v, d, d)
    # muMat <- matrix(mu, d, d)
    # vNorm <- sqrt(sum(v^2) / 2)
    # res <- diag(3) + sin(vNorm) / vNorm * vMat + (1 - cos(vNorm)) / vNorm^2 * vMat %*% vMat
    # res <- as.numeric(res %*% muMat)
  # }

  # res
# }



# Get the nearest orthogonal matrix
NearestOrth <- function(A) {
  if (length(A) == 1) {
    A <- matrix(A)
  }
  stopifnot(is.matrix(A))
  stopifnot(nrow(A) == ncol(A))

  s <- svd(A)

  tcrossprod(s[['u']], s[['v']])
}


# Get the nearest SO matrix
#' @export
project.SO <- function(mfd, A) {

  A <- as.matrix(A)
  d <- round(sqrt(nrow(A)))

  res <- apply(A, 2, function(a) {
    a <- matrix(a, d, d)
    N <- NearestOrth(a)
    d <- nrow(N)
    if (!isSO(N)) {
      N[, d] <- -N[, d]
    }
    c(N)
  })
  res
}


# Calculate the parameter n of SO(n) given the length of the lower triangle
calcGeomPar.SO <- function(mfd, dimIntrinsic, dimAmbient, dimTangent) {

  if (!missing(dimIntrinsic)) {
    dimTangent <- calcTanDim.SO(dimIntrinsic=dimIntrinsic)
  } else if (!missing(dimAmbient)) {
    dimTangent <- calcTanDim.SO(dimAmbient=dimAmbient)
  } else if (!missing(dimTangent)) {
  }

  d0 <- sqrt(2 * dimTangent + 1/4) + 1/2 # Root of a quadratic equation
  n <- round(d0)

  if (n != d0) {
    stop('The input is not the dimTangent of a lower trianglular matrix')
  }

  n
}


calcIntDim.SO <- function(mfd, dimAmbient, dimTangent) {

  if (!missing(dimAmbient)) {
    n <- round(sqrt(dimAmbient))
  } else if (!missing(dimTangent)) {
    n <- calcGeomPar.SO(dimTangent=dimTangent)
  }
  n * (n - 1) / 2
 
}


calcTanDim.SO <- function(mfd, dimAmbient, dimIntrinsic) {

  if (!missing(dimAmbient)) {
    calcIntDim.SO(dimAmbient=dimAmbient)
  } else if (!missing(dimIntrinsic)) {
    dimIntrinsic
  }

}


MakePhi.SO <- 
  function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           type=c('cos', 'sin', 'fourier', 'legendre01', 'poly')) {

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


NoiseTangent.SO <- function(mfd, n, mu, sigma2=0, epsFun = rnorm) {

  m <- ncol(mu)
  dimAmbient <- nrow(mu)
  dimTangent <- calcTanDim.SO(mfd, dimAmbient)
  res <- array(epsFun(n * dimTangent * m), c(n, dimTangent, m)) * sqrt(sigma2)
  res

}
