# devtools::load_all()
library(testthat)

test_that('MakeMfdProcess and RFPCA works for L2 manifold', {

  n <- 50
  m <- 21
  p <- 40
  tpts <- seq(0, 1, length.out=m)
  spts <- seq(0, 1, length.out=p)
  mfdL2 <- structure(1, class='L2')
  sparsity <- p
  K <- 6
  lambda <- 0.05 ^ (seq_len(K) / 3)

  # mu is a beta pdf 
  a <- seq(1, 3, length.out=m)
  b <- 1 + 2 * tpts * (1 - tpts)
  mu <- vapply(seq_len(m), function(i) {
    res <- dbeta(spts, a[i], b[i]) + 1
    sqrt(res / sum(res) * p)
  }, rep(0, p))

  set.seed(1)
  samp <- MakeMfdProcess(mfdL2, n, mu, tpts, K = K, lambda=lambda)

  expect_equal(apply(samp$phi, 2, function(x) mean(x^2)), rep(1, K), tolerance=1e-1)
  expect_equal(crossprod(samp$phi) / m / p, diag(K), tolerance=1e-1)
  # matplot(t(samp$X[3, , ]), type='l')
  # matplot(t(samp$XNoisy[3, , ]), type='l')
  # image(samp$X[1, , ])
  spSamp <- SparsifyM(samp$X, samp$T, m)
  yList <- spSamp$Ly
  tList <- spSamp$Lt

  resL2 <- RFPCA(yList, tList, list(maxK=K, mfd=mfdL2, dataType='Dense'))

  phiMat <- matrix(resL2$phiObsTrunc, m * p, K)
  expect_equal(crossprod(phiMat) / (m - 1) / (p - 1), diag(K), tolerance=1e-1)
  expect_equal(abs(crossprod(phiMat[, 1:5], samp$phi[, 1:5]) / m / p), diag(5), tolerance=2e-1)
  expect_equal(resL2$lam, lambda, tolerance=1e-1)
  expect_equal(abs(diag(cov(resL2$xi, samp$xi))), lambda, tolerance=1e-1)
  expect_true(all(abs(diag(cor(resL2$xi, samp$xi)))[1:5] > 0.8))
})


test_that('In the noiseless situation, using all K results in near perfect fitted trajectories', {
  n <- 3
  # m <- 51
  for (m in c(1, 51)) {
    p <- 4
    tpts <- seq(0, 1, length.out=m)
    spts <- seq(0, 1, length.out=p)
    mfdL2 <- structure(1, class='L2')
    sparsity <- p
    K <- 3
    lambda <- 0.05 ^ (seq_len(K) / 3)

    mu <- matrix(0, p, m)

    set.seed(1)
    samp <- MakeMfdProcess(mfdL2, n, mu, tpts, K = K, lambda=lambda)
    samp$X <- abind::abind(samp$X, -samp$X, along=1)
    samp$xi <- abind::abind(samp$xi, -samp$xi, along=1)

    spSamp <- SparsifyM(samp$X, samp$T, m)
    yList <- spSamp$Ly
    tList <- spSamp$Lt

    resL2 <- RFPCA(yList, tList, list(mfd=mfdL2, dataType='Dense', maxK=K))
    errFitL2 <- ErrFitted(mfdL2, samp$X, resL2)

    expect_equal(errFitL2[K], 0)
  }
})

test_that('MakeMfdProcess and RFPCA works for DensL2 manifold', {

  n <- 50
  nnew <- 100
  for (m in c(1, 21)) {
    p <- 40
    tpts <- seq(0, 1, length.out=m)
    spts <- seq(0, 1, length.out=p)
    mfdHS <- structure(1, class=c('HS', 'L2')) # Square-root densities
    mfdDens <- structure(1, class=c('Dens', 'L2')) # Densities, as elements in L2
    sparsity <- p
    K <- 15
    lambda <- 0.05 ^ (seq_len(K) / 3)

    # mu is a beta pdf 
    a <- seq(1.2, 3, length.out=m)
    b <- 1 + 2 * tpts * (1 - tpts)
    mu <- vapply(seq_len(m), function(i) {
      res <- dbeta(spts, a[i], b[i]) + 1
      sqrt(res / sum(res) * (length(res) - 1))
    }, rep(0, p))

    # The numerical precision for the left endpoint rule for the legendre basis is poor. Don't test at this time
    for (bt in c('cos', 'sin', 'fourier')) {
      phi <- MakePhi.HS(mfdHS, K, mu, tpts, bt)
      phiHS <- MakePhi.HS(mfdHS, K, mu, tpts, bt)
      # phi1 <- MakePhi.HS(mfdHS, K, mu, tpts, bt)
      phiMat <- matrix(phi, ncol=K)
      # phi are orthonormal
      expect_equal(crossprod(phiMat) / p / m, diag(K), tolerance=1e-1)

      # phi(t) are orthogonal to mu(t)
      tmp <- apply(phi, 3, function(phi1) {
        colSums(t(phi1) * mu) / p / m
      })
      tmp <- matrix(tmp, m, K)
      expect_equal(unname(tmp), matrix(0, m, K))
    }

    set.seed(1)
    samp <- MakeMfdProcess(mfdHS, n, mu, tpts, K = K, lambda=lambda)
    sampNew <- MakeMfdProcess(mfdHS, nnew, mu, tpts, K = K, lambda=lambda)
    # Samples are on HS
    expect_equal(unname(apply(samp$XNoisy, c(1, 3), function(x) sum(x^2) / (p - 1))), matrix(1, n, m))
    expect_equal(unname(apply(samp$X, c(1, 3), function(x) sum(x^2) / (p - 1))), matrix(1, n, m))
    # Variance of xi are near lambda
    expect_equal(diag(var(samp$xi)), lambda[seq_len(K)], tolerance=1e-1)

    phiAr <- array(samp$phi, c(m, p, K))
    # phi_j(t) are orthogonal to mu(t)
    expect_equal(matrix(apply(phiAr, 3, function(phi) colMeans(t(phi) * mu)), m, K), 
                 matrix(0, m, K))
    # phi_j are orthogonal
    expect_equal(crossprod(samp$phi) / m / p, diag(K), tolerance=1e-1)

    # image(samp$X[7, , ])
    # image(phiAr[, , 3])
    spSamp <- SparsifyM(samp$X, samp$T, m)
    yList <- spSamp$Ly
    y2List <- sapply(yList, `^`, 2, simplify=FALSE)
    tList <- spSamp$Lt
    spSampNew <- SparsifyM(sampNew$X, sampNew$T, m)
    yListNew <- spSampNew$Ly
    tListNew <- spSampNew$Lt

    resHS <- RFPCA(yList, tList, list(maxK=K, mfd=mfdHS, dataType='Dense'))
    resDensL2 <- RFPCA(y2List, tList, list(maxK=K, mfd=mfdDens, dataType='Dense'))

    fitObs <- fitted(resHS, grid='obs')
    fitWork <- fitted(resHS, grid='work')
    predObs <- predict(resHS, yList, tList, type='traj', xiMethod='IN')
    xiPred <- predict(resHS, yList, tList, type='xi', xiMethod='IN')
    expect_equal(dim(fitObs)[3], m)
    expect_equal(dim(fitWork)[3], length(resHS$workGrid))
    expect_equal(unname(fitObs), predObs)
    expect_equal(unname(resHS$xi), xiPred)

    expect_equal(resHS$lam, lambda[seq_len(K)], tolerance=1e-1)
    expect_equal(resDensL2$lam, lambda[seq_len(K)], tolerance=1e-1)
    expect_equal(abs(diag(cov(resHS$xi, samp$xi))), lambda, tolerance=1e-1)
    expect_true(all(abs(diag(cor(resHS$xi, samp$xi)))[1:5] > 0.8))

    errFitHS <- ErrFitted(mfdHS, samp$X, resHS)
    errPredHS <- ErrFitted(mfdHS, sampNew$X, resHS, yListNew, tListNew)
    errPredExt <- ErrFitted(mfdHS, sampNew$X, resHS, yListNew, tListNew, expMap=FALSE)
    expect_true(all(errPredHS > errFitHS - 1e-15))
    expect_true(all(errPredExt > errPredHS - 1e-15))
  }
})
