# devtools::load_all()
library(testthat)

# Some simulation
test_that('RFPCA works for sparse case', {
  n <- 50
  m <- 20
  K <- 20
  KDirect <- 4
  lambda <- 0.07 ^ (seq_len(K) / 2)
  p <- 3
  basisType <- 'legendre01'
  sparsity <- 10:15
  sigma2 <- 0.01
  if (p == 3) {
    muList <- list(
      function(x) x * 2, 
      function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
      function(x) rep(0, length(x))
    )
  } else if (p == 4) {
    muList <- list(
      function(x) x * sqrt(2), 
      function(x) x * sqrt(2), 
      function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
      function(x) rep(0, length(x))
    )
  }
  pts <- seq(0, 1, length.out=m)
  mfd <- structure(1, class='Sphere')
  mu <- Makemu(mfd, muList, c(rep(0, p - 1), 1), pts)
  expect_equal(colSums(mu^2), rep(1, ncol(mu)))

  # Generate noiseless samples
  # CreateBasis <- fdapace::CreateBasis
  set.seed(1)
  samp <- MakeSphericalProcess(n, mu, pts, K = K, lambda=lambda, basisType=basisType, sigma2=sigma2)
  spSamp <- SparsifyM(samp$X, samp$T, sparsity)
  yList <- spSamp$Ly
  tList <- spSamp$Lt

  # smooth
  bw <- 0.2
  kern <- 'epan'
  assumeError <- TRUE

  resSp <- RFPCA(yList, tList, list(userBwMu=bw, userBwCov=bw * 2, kernel=kern, maxK=K))
  resEu <- RFPCA(yList, tList, list(userBwMu=bw, userBwCov=bw * 2, kernel=kern, maxK=K, mfdName='euclidean'))
  KCFPCA <- 10
  resC <- CFPCA(yList, tList, list(userBwMu=bw, userBwCov=bw * 2, kernel=kern, KUse=KCFPCA, FVEthreshold=1))

  expect_equal(c(resSp$muObs), c(mu), tolerance=0.2)
  expect_equal(c(resEu$muObs), c(mu), tolerance=0.2)
  expect_equal(resSp$lam, lambda, tolerance=0.02)
  expect_equal(resEu$lam, lambda, tolerance=0.02)
  expect_equal(resC$lam[seq_len(KCFPCA)], lambda[seq_len(KCFPCA)], tolerance=0.1)
  expect_equal(resSp$sigma2, sigma2, tolerance=0.03)
  expect_equal(resEu$sigma2, sigma2, tolerance=0.03)
})


test_that('RFPCA works for dense case', {
  n <- 100
  m <- 50
  K <- 20
  KDirect <- 4
  lambda <- 0.07 ^ (seq_len(K) / 2)
  dimTangent <- 80
  basisType <- 'legendre01'
  sparsity <- m
  sigma2 <- 0
  muScale <- 1
  muList <- c(
    rep(list(function(x) x * 2 * muScale / sqrt(dimTangent - 2)), dimTangent - 2),
    function(x) sin(x * 1 * pi) * pi / 2 * 0.6 * muScale,
    function(x) rep(0, length(x))
  )
  pts <- seq(0, 1, length.out=m)
  mfd <- structure(1, class='Sphere')
  mu <- Makemu(mfd, muList, c(rep(0, dimTangent - 1), 1), pts)
  expect_equal(colSums(mu^2), rep(1, ncol(mu)))

  # Generate noiseless samples
  # CreateBasis <- fdapace::CreateBasis
  set.seed(1)
  samp <- MakeSphericalProcess(n, mu, pts, K = K, lambda=lambda, basisType=basisType, sigma2=sigma2)
  spSamp <- SparsifyM(samp$X, samp$T, sparsity)
  yList <- spSamp$Ly
  tList <- spSamp$Lt

  # print(system.time({
  # print(profvis::profvis(
    resSp <- RFPCA(yList, tList, list(maxK=K, mfd=mfd))
  # ))
  # }))
  # print(system.time({
    resEu <- RFPCA(yList, tList, list(maxK=K, mfdName='euclidean'))
  # }))

  expect_equal(c(resSp$muObs), c(mu), tolerance=0.1)
  expect_equal(c(resEu$muObs), c(mu), tolerance=0.1)
  expect_equal(resSp$lam, lambda, tolerance=0.02)
  expect_equal(resEu$lam, lambda, tolerance=0.02)
})
