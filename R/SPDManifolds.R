## Classes for general symmetric positive definite matrix manifolds, with different metrics (log-Euclidean & affinite-invariant). c.f. Arsigny et al 2007
# The tangent space is the space of symmetric matrices. This space is represented by a p-by-p symmetric matrix vectorized. Each basis vector has all 0s but either 1. has a single 1 entry in the diagonal, or 2. two 1/sqrt{2} in a pair of off-diagonal symmetric entries. 
#
# Log-Euclidean:
# Choose basis (chart) such that matrices are represented by their matrix-logs, and the metric is given by the Euclidean inner product
#
# Affine-invariant
# Choose basis such that the logarithm map log_{S_1}(S_2) is represented by 
# L'=log(S_1^{-1/2} S_2 S_1^{-1/2})
# The exponential map exp_{S_1}(L') is then
# S_2 = S_1^{1/2}exp(L')S_1^{1/2}


#' @export
#' @describeIn metric Method
metric.LogEu <- function(mfd, p, U, V) {
  metric.default(mfd, p, U, V)
}


#' @export
#' @describeIn norm Method
norm.LogEu <- function(mfd, p, U) {
  norm.default(mfd, p, U)
}


# #' Geodesic distance of points on the manifold
# #' 
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param X A D*n matrix. Each column represents a point on manifold
# #' @param Y A D*n matrix. Each column represents a point on manifold
# #' @return A 1*n vector. The \emph{i}th element is the geodesic distance of \code{X[, i]} and \code{Y[, i]}

#' @export
#' @describeIn distance Method
distance.LogEu <- function(mfd, X, Y, assumeLogRep=FALSE) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  stopifnot(nrow(X) == nrow(Y))
  if (length(X) == 0 || length(Y) == 0) {
    return(numeric(0))
  }

  if (assumeLogRep) {
    return(distance.Euclidean(X=X, Y=Y))
  } else if (!assumeLogRep) {
    d <- round(sqrt(nrow(X)))

    if (ncol(X) == 1) {
      logX <- matrix(LogM(matrix(X, d, d)), d^2, ncol(Y))
    } else {
      logX <- matrix(apply(X, 2, function(x) {
        LogM(matrix(x, d, d))
      }), d^2, ncol(X))
    }

    if (ncol(Y) == 1) {
      logY <- matrix(LogM(matrix(Y, d, d)), d^2, ncol(X))
    } else {
      logY <- matrix(apply(Y, 2, function(x) {
        LogM(matrix(x, d, d))
      }), d^2, ncol(Y))
    }

    return(norm.default(U=logX - logY))
  }
}


# #' Riemannian exponential map at a point

#' @export
#' @describeIn rieExp Method
rieExp.LogEu <- function(mfd, p, V, tol=1e-10) {

  V <- as.matrix(V)

  if (length(V) == 0 || (!missing(p) && length(p) == 0)) {
    return(matrix(0, nrow(V), 0))
  }

  n <- ncol(V)
  d <- calcGeomPar.SPD(dimTangent=nrow(V))


  if (missing(p)) {
    logp <- matrix(0, d^2, n)
  } else {
    if (is.matrix(p) && ncol(p) == n) {
      logp <- apply(p, 2, logmvec, d=d)
    } else {
      logp <- matrix(logmvec(p, d), length(p), n)
    }
  }

  res <- apply(logp + V, 2, expmvec, d=d)
  matrix(res, d^2, n)
}


#' @export
#' @describeIn rieLog Method
rieLog.LogEu <- function(mfd, p, X, tol=1e-10) {

  X <- as.matrix(X)
  d2 <- nrow(X)
  n <- ncol(X)
  d <- round(sqrt(d2))

  if (missing(p)) {
    p <- matrix(diag(nrow=d))
    logp <- matrix(0, d^2, n)
  } else {
    p <- as.matrix(p)
  }
  stopifnot(nrow(p) == nrow(X))

  if (length(p) == 0 || length(X) == 0) {
    return(matrix(0, nrow(p), 0))
  }

  if (ncol(p) == 1) {
    n <- ncol(X)
    logp <- matrix(logmvec(p, d), d2, n)
    logX <- matrix(apply(X, 2, logmvec, d=d), d2, n)
  } else if (ncol(X) == 1) {
    n <- ncol(p)
    logX <- matrix(logmvec(X, d), d2, n)
    logp <- matrix(apply(p, 2, logmvec, d=d), d2, n)
  } else if (ncol(X) == ncol(p)) {
    n <- ncol(p)
    logp <- matrix(apply(p, 2, logmvec, d=d), d2, n)
    logX <- matrix(apply(X, 2, logmvec, d=d), d2, n)
  }

  logX - logp
}


#' @export
#' @describeIn metric Method
metric.AffInv <- function(mfd, p, U, V) {
  metric.default(mfd, p, U, V)
}


#' @export
#' @describeIn norm Method
norm.AffInv <- function(mfd, p, U) {
  norm.default(mfd, p, U)
}


# A vector X representing a single point would be faster than using a matrix X due to the rieLog.AffInv

#' @export
#' @describeIn distance Method
distance.AffInv <- function(mfd, X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  stopifnot(nrow(X) == nrow(Y))
  distAffInv(X, Y) # Requires individual matrices combined into columns
}



#' @export
#' @describeIn rieExp Method
rieExp.AffInv <- function(mfd, p, V, tol=1e-10) {

  if (length(p) == 0 || length(V) == 0) {
    if (is.matrix(p)) {
      return(matrix(0, nrow(V), 0))
    } else {
      return(matrix(0, length(p), 0))
    }
  }

  p <- as.matrix(p)
  V <- as.matrix(V)
  stopifnot(nrow(p) == nrow(V))

  n <- ncol(V)
  d <- calcGeomPar.SPD(dimTangent=nrow(V))

  if (!is.matrix(p) || ncol(p) == 1) {
    # Common base point
    pHalf <- SqrtM(matrix(p, d, d), 1/2)
    res <- apply(V, 2, function(v) {
      affineCenter(s2=c(ExpM(matrix(v, d, d))), S1HalfInv=pHalf, d=d)
    })
  } else if (is.matrix(p) && ncol(p) == n) { 
    # Each tangent vector has a different base point
    res <- vapply(seq_len(ncol(V)), function(i) {
      affineCenter(p[, i], c(ExpM(matrix(V[, i], d, d))), d=d, power=1/2)
    }, rep(0, d^2))
  }

  matrix(res, d^2, n)
}



#' @export
#' @describeIn rieLog Method
rieLog.AffInv <- function(mfd, p, X, tol=1e-10) {

  p <- as.matrix(p)
  X <- as.matrix(X)
  stopifnot(nrow(p) == nrow(X))
  logAffInv(p, X)

}


# Take log of a matrix.
# x Vectorized matrix
# d Dimension of the matrix
logmvec <- function(x, d) {
  LogM(matrix(x, d, d))
}


expmvec <- function(x, d) {
  ExpM(matrix(x, d, d))
}


# Square-root matrix, through eigendecomposition. Assume is symmetric
SqrtM <- function(mat, power=1/2, tol=100 * .Machine$double.eps) {

  mat <- as.matrix(mat)

  # if (!isSymmetric(mat, tol)) {
      # stop('`mat` should be symmetric')
  # }

  mat <- as.matrix(mat)
  eigRes <- eigen(mat, symmetric=TRUE)
  # eigRes[['values']][eigRes[['values']] < 0] <- 0

  newEigVal <- eigRes[['values']]^power

  res <- eigRes[['vectors']] %*% 
    diag(newEigVal, length(newEigVal)) %*% 
    t(eigRes[['vectors']])
  res
}


# s1, s2 are vectors
# returns a d by d symmetric matrix 
affineCenter <- function(s1, s2, S1HalfInv, d, power=-1/2) {
  S2 <- matrix(s2, d, d)
  if (missing(S1HalfInv)) {
    S1 <- matrix(s1, d, d)
    S1HalfInv <- SqrtM(S1, power=power)
  }
  MakeSym(S1HalfInv %*% S2 %*% t(S1HalfInv))
}


# n: Dimension of a matrix
# cf https://math.stackexchange.com/questions/1143614/is-matrix-transpose-a-linear-transformation
TransposeAsMatrix <- function(n) {
  ll <- lapply(seq_len(n), function(i) {
    a <- rep(0, n)
    a[i] <- 1
    l <- rep(list(matrix(a, 1, n)), n)
    Matrix::bdiag(l)
  })
  do.call(rbind, ll)
}


## Helper functions for SPD matrices
isSPSD <- function(x, tol=1e-14) {
  if (length(x) == 1) {
    x <- matrix(x, 1, 1)
  }

  ev <- eigen(x)$values
  isSymmetric(x, tol) && min(ev) >= 0
}


isSPD <- function(x, tol=1e-14) {
  if (length(x) == 1) {
    x <- matrix(x, 1, 1)
  }

  ev <- eigen(x)$values
  isSymmetric(x, tol) && min(ev) >= tol
}


#' Make a symmetric matrix by specifying a near-symmetric matrix M, or the lower trianglular elements lowerTri with diagonal.
#' @param M A near-symmetric matrix
#' @param lowerTri A vector containing the lower triangular elements of the matrix. This is an alternative way to specify the matrix. 
#' @param doubleDiag Only meaningful for lowerTri is not missing. Whether the diagonal elements should be multiplied by sqrt(2) for doubling the squared norm of the lower triangle
#' @export
MakeSym <- function(M, lowerTri, doubleDiag=FALSE) {

  if (!missing(M)) {
    A <- (M + t(M)) / 2
  } else if (!missing(lowerTri)) {
    p <- length(lowerTri)
    d <- calcGeomPar.SPD(dimIntrinsic=p)

    A <- matrix(0, d, d)
    A[lower.tri(A, diag=TRUE)] <- lowerTri
    A <- A + t(A) - diag(diag(A), d)
    if (doubleDiag) {
      diag(A) <- diag(A) * sqrt(2)
    } 
  }

  A
}

# Get the nearest positive semidefinite matrix
NearestSPSD <- function(A) {
  n <- sqrt(length(A))
  if (n == 1) {
    A <- matrix(A)
  }
  stopifnot(is.matrix(A))
  stopifnot(nrow(A) == ncol(A))

  eig <- eigen(A)
  eig$values[eig$values < 0] <- 0

  eig$vectors %*% diag(eig$values, n) %*% t(eig$vectors)
}


#' @export
#' @describeIn project Method
project.AffInv <- function(mfd, A) {
  NextMethod('project')
}


#' @export
#' @describeIn project Method
project.LogEu <- function(mfd, A) {
  NextMethod('project')
}


#' @export
#' @describeIn project Method
# Get the nearest positive semidefinite matrix
project.SPD <- function(mfd, A) {

  A <- as.matrix(A)
  d <- round(sqrt(nrow(A)))

  res <- apply(A, 2, function(a) {
    a <- matrix(a, d, d)
    N <- NearestSPSD(a)
    c(N)
  })
  matrix(res, ncol=ncol(A))
}


# Project extrinsic tangent space data onto the tangent space at p
# p must be a single point on the manifold
#' @export
#' @describeIn projectTangent Method
projectTangent.SPD <- function(mfd, p, X, projMatOnly=FALSE) {

  if (!missing(p)) {
    n <- length(p)
    p <- as.matrix(p)
  } else if (!missing(X)) {
    X <- as.matrix(X)
    n <- nrow(X)
  } else {
    stop('p and X cannot both be missing')
  }
  d <- round(sqrt(n))

  if (projMatOnly) {
    projMat <- (diag(nrow=d ^ 2) + TransposeAsMatrix(d)) / 2
    return(projMat)
  } else {
    res <- apply(X, 2, function(x) {
      xM <- matrix(x, d, d)
      (xM + t(xM)) / 2
    })
    matrix(res, d ^ 2, ncol(X))
  }
}


# metric.SPD <- function(mfd, p, U, V) {
  # NextMethod('metric')
# }


# norm.SPD <- function(mfd, p, U) {
  # NextMethod('norm')
# }


# Calculate the parameter n of SPD(n) given the length of the lower triangle with diagonal. 
# Let l=(n + 1) * n / 2

#' @export
calcGeomPar.SPD <- function(mfd, dimIntrinsic, dimAmbient, dimTangent) {

  if (!missing(dimAmbient)) {
    d0 <- sqrt(dimAmbient)
  } else if (!missing(dimTangent)) {
      d0 <- sqrt(dimTangent)
  } else if (!missing(dimIntrinsic)) {
    l <- dimIntrinsic
    d0 <- sqrt(2 * l + 1/4) - 1/2 # Root of a quadratic equation
  }

  n <- round(d0)

  if (n != d0) {
    stop('The input does not correspond to a SPD matrix')
  }

  as.integer(n)
}


#' @export
calcIntDim.SPD <- function(mfd, geomPar, dimAmbient, dimTangent) {

  if (!missing(geomPar)) {
    return(geomPar * (geomPar + 1) / 2)
  } else if (!missing(dimAmbient)) {
    dd <- dimAmbient
  } else if (!missing(dimTangent)) {
    dd <- dimTangent
  }
  d0 <- sqrt(dd)
  n <- round(d0)
  stopifnot(d0 == n)
  as.integer(n * (n + 1) / 2)
}


#' @export
calcTanDim.SPD <- function(mfd, geomPar, dimAmbient, dimIntrinsic) {

  if (!missing(geomPar)) {
    as.integer(geomPar^2)
  } else if (!missing(dimAmbient)) {
    as.integer(dimAmbient)
  } else if (!missing(dimIntrinsic)) {
    n <- calcGeomPar.SPD(mfd, dimIntrinsic=dimIntrinsic)
    as.integer(n^2)
  }

}


#' @export
calcAmbDim.SPD <- function(mfd, geomPar, dimIntrinsic, dimTangent) {
  if (!missing(geomPar)) {
    as.integer(geomPar ^ 2)
  } else if (!missing(dimIntrinsic)) {
    as.integer(calcGeomPar.SPD(mfd, dimIntrinsic=dimIntrinsic)^2)
  } else if (!missing(dimTangent)) {
    as.integer(calcGeomPar.SPD(mfd, dimTangent=dimTangent)^2)
  }
}


MakePhi.SPD <- 
  function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           type=c('cos', 'sin', 'fourier', 'legendre01', 'poly')) {

  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }
  dimAmbient <- nrow(mu)
  dimIntrinsic <- calcIntDim.SPD(mfd, dimAmbient = dimAmbient)
  dimTangent <- calcTanDim.SPD(mfd, dimAmbient = dimAmbient)
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

  ptsSqueeze <- do.call(c, lapply(seq_len(dimIntrinsic), function(i) {
                          pts + (max(pts) + mean(diff(pts))) * (i - 1)
                        }))
  ptsSqueeze <- ptsSqueeze / max(ptsSqueeze)
  
  phi <- array(fdapace:::CreateBasis(K, ptsSqueeze, type), c(m, dimIntrinsic, K)) / 
    sqrt(dimIntrinsic)
  # Generate the lower triangle, and then symmetrize
  tmp <- apply(phi, c(1, 3), function(x) {
    x <- MakeSym(lowerTri = x, doubleDiag=TRUE)
    x 
  }) / sqrt(2)
  phi <- aperm(tmp, c(2, 1, 3))

  dimnames(phi) <- list(t=pts, j=seq_len(dimTangent), k=seq_len(K))
  phi
}


#' @export
#' @describeIn NoiseTangent Method
NoiseTangent.SPD <- function(mfd, n, mu, sigma2=0, epsFun = rnorm) {

  m <- ncol(mu)
  dimAmbient <- nrow(mu)
  dimIntrinsic <- calcIntDim.SPD(mfd, dimAmbient=dimAmbient)
  res <- array(epsFun(n * dimIntrinsic * m), c(n, dimIntrinsic, m)) * sqrt(sigma2 / 2)
  tmp <- apply(res, c(1, 3), function(x) MakeSym(lowerTri=x, doubleDiag=TRUE))
  res <- aperm(tmp, c(2, 1, 3))

  res

}


# Use a trick to convert a matrix representation to the log-lower trianglar representation
# l A matrix whose columns are vectorized matrices
ToLower <- function(l, takeLog=TRUE) {
  l <- as.matrix(l)
  d <- round(sqrt(nrow(l)))
  a <- apply(l, 2, function(x) {
               if (takeLog) {
                 m <- logmvec(x, d=d)
               } else {
                 m <- matrix(x, d, d)
               }
               diag(m) <- diag(m) / sqrt(2)
               m[lower.tri(m, diag=TRUE)]
  })

  a
}


#' @export
#' @describeIn origin The origin is the identity matrix.
origin.SPD <- function(mfd, dimIntrinsic) {
  # browser()
  as.matrix(c(diag(calcGeomPar.SPD(mfd, dimIntrinsic=dimIntrinsic))))
}


#' @export
#' @describeIn basisTan Method
basisTan.SPD <- function(mfd, p) {

  # browser()
  d <- round(sqrt(length(p)))
  dInt <- RFPCA::calcIntDim(mfd, geomPar=d)
  P <- matrix(p, d, d)
  ind <- which(lower.tri(P, diag=TRUE), arr.ind=TRUE, useNames=FALSE)
  res <- apply(ind, 1, function(ii) {
    M <- matrix(0, d, d)
    if (ii[1] == ii[2]) {
      M[ii[1], ii[2]] <- 1
    } else {
      M[ii[1], ii[2]] <- 1/sqrt(2)
      M[ii[2], ii[1]] <- 1/sqrt(2)
    }
    c(M)
  })
  res <- matrix(res, nrow=d^2)
  res
}


# basisTan.AffInv <- function(mfd, p) {
  # NextMethod()
# }


# basisTan.LogEu <- function(mfd, p) {
  # NextMethod()
# }


