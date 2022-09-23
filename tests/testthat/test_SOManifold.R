library(testthat)

test_that('MakeMfdProcess, GCVFrechetMeanCurve, frechetMeanCurve works for SO', {

  n <- 200
  p <- 3
  m <- 50
  pts <- seq(0, 1, length.out=m)
  sparsity <- 2:5
  bwMu <- 0.1
  kern <- 'epan'
  muScale <- 1

  if (p == 1) {
    muList <- list( function(x) x * muScale)
  } else if (p == 3) {
    muList <- list(
      function(x) x * 2 * muScale, 
      function(x) sin(x * 1 * pi) * pi / 2 * 0.6 * muScale,
      function(x) rep(0, length(x))
    )
  } else if (p == 4) {
    # stop()
    # muList <- list(
      # function(x) x * sqrt(2), 
      # function(x) x * sqrt(2), 
      # function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
      # function(x) rep(0, length(x))
    # )
  }

  set.seed(1)
  mfd <- structure(list(), class='SO')
  d <- calcGeomPar(mfd, dimTangent = p)
  mu <- Makemu(mfd, muList, as.numeric(diag(d)), pts)
  samp <- MakeMfdProcess(mfd, n, mu, pts, K = 3, lambda=c(0.3, 0.2, 0.1), xiFun=rnorm, epsFun=rnorm, sigma2=0.01)
  spSamp <- SparsifyM(samp$XNoisy, samp$T, sparsity)

  yList <- spSamp$Ly
  tList <- spSamp$Lt
  Ymat <- do.call(cbind, yList)
  Tvec <- do.call(c, tList)
  uniqT <- sort(unique(Tvec))
  tInd <- match(uniqT, pts)
  ord <- order(Tvec)

  # a <- profvis({
  muM <- frechetMeanCurve(mfd, bwMu, kern, xin=Tvec[ord], yin=Ymat[, ord], xout=pts, npoly=0)
  # })
  bwMuGCV <- GCVFrechetMeanCurve(mfd, yList, tList, kern, 0, 'Sparse')$bOpt
  bwMuGCVEu <- GCVFrechetMeanCurve(structure(1, class='Euclidean'), yList, tList, kern, 0, 'Sparse')$bOpt
  matplot(t(mu), type='l', lty=1)
  matplot(t(muM), type='l', lty=2)
  expect_equal(mu, muM, tolerance=2e-1)
  expect_equal(bwMuGCV, 0.1, tolerance=0.1)
  expect_equal(bwMuGCVEu, 0.1, tolerance=0.1)

})


test_that('Makemu and Makephi works for SO and sphere', {

  mfd <- structure(list(), class='SO')
  m <- 50
  pts <- seq(0, 1, length.out=m)

  v0 <- function(x) rep(0, length(x))
  v1 <- function(x) x * 2 
  v2 <- function(x) exp(- (x * 2 - 1/4)^2 * 3) - exp(- 3 / 16)
  v3 <- function(x) rep(0, length(x))

  res0 <- Makemu(mfd, list(v0), diag(2), pts)
  expect_equal(res0, matrix(c(1, 0, 0, 1), 4, m))

  res1 <- Makemu(mfd, list(v1), diag(2), pts)
  res2 <- Makemu(mfd, list(v0, v1, v2), diag(3), pts)

  # The first row is the identity
  expect_equal(res1[, 1], as.numeric(diag(2)))
  expect_equal(res2[, 1], as.numeric(diag(3)))

  # The mean function lies on the manifold
  expect_equal(t(apply(res1, 2, function(x) crossprod(matrix(x, 2, 2)))), 
               matrix(diag(2), m, 2^2, byrow=TRUE))
  expect_equal(t(apply(res2, 2, function(x) crossprod(matrix(x, 3, 3)))), 
               matrix(diag(3), m, 3^2, byrow=TRUE))

  mfdS <- structure(list, class='Sphere')
  expect_equal(colSums(Makemu(mfdS, list(v1, v2, v3), c(0, 0, 1))^2),
               rep(1, m))

  d <- 2
  phi1 <- MakephiSO(1, t(res1), type='sin')
  phi2 <- MakephiSO(3, t(res1))
  phi3 <- MakephiSO(3, t(res2))

  # phi_k(t) is skew symmetric
  expect_true(all(apply(phi1, c(1, 3), function(v) isSkewSym(matrix(v, 2, 2)))))
  expect_true(all(apply(phi2, c(1, 3), function(v) isSkewSym(matrix(v, 2, 2)))))
  expect_true(all(apply(phi3, c(1, 3), function(v) isSkewSym(matrix(v, 3, 3)))))
  # phi_k(t) integrates to 1
  expect_equal(max(abs(apply(phi1, 3, function(phi) sum(phi^2)) / m - 1)), 0, tolerance=5e-2)
  expect_equal(max(abs(apply(phi2, 3, function(phi) sum(phi^2)) / m - 1)), 0, tolerance=5e-2)
  expect_equal(max(abs(apply(phi3, 3, function(phi) sum(phi^2)) / m - 1)), 0, tolerance=5e-2)


})

