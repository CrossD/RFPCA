## Generate process on a sphere by exponential maps
# Return a three dimensional array X[i, j, t] = X_{ij}(t)
# sigma2: The entriwise (isotropic) variance for added errors
# xiFun: A function to generate xi. The first argument is the number of samples to generate.
# epsFun: A function to generate noise, similar to xiFun.
# epsBase: What is the base point p for adding noises. The noisy observations are Exp_p(eps), where p = 'mu' or 'X'
MakeSphericalProcess <- function(
  n, mu, pts=seq(0, 1, length.out=ncol(mu)), K=2, lambda=rep(1, K), sigma2=0, 
  xiFun = stats::rnorm, epsFun = stats::rnorm, epsBase = c('mu', 'X'), 
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



