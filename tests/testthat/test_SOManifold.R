devtools::load_all()
library(testthat)

test_that('Helper functions for SO(n) works', {
  
  # isSkewSym
  expect_true(isSkewSym(0))
  expect_false(isSkewSym(1))
  expect_true(isSkewSym(matrix(c(0, 1, -1, 0), 2, 2)))
  expect_false(isSkewSym(matrix(c(0, 1, 1, 0), 2, 2)))
  expect_true(isSkewSym(0.01, tol=0.1))
  expect_false(isSkewSym(0.1, tol=0.01))

  # isOrth
  expect_true(isOrth(1))
  expect_false(isOrth(0))
  expect_true(isOrth(diag(2)))
  expect_true(isOrth(diag(c(1, -1))))
  expect_false(isOrth(diag(2) + 0.1))

  # isSO
  expect_true(isSO(diag(2)))
  expect_false(isSO(diag(c(1, -1))))

  # NearestOrth
  set.seed(1)
  d <- 3
  A <- matrix(rnorm(d^2), d, d)
  O1 <- NearestOrth(A)
  expect_true(isOrth(O1))

  # NearestSO
  set.seed(1)
  d <- 3
  A <- matrix(rnorm(d^2), d, d)
  B <- diag(d)
  SO1 <- project.SO(A=cbind(c(A), c(B)))
  expect_true(isSO(matrix(SO1[, 1], d, d)))
  expect_equal(SO1[, 2], c(B))

  # MakeSkewSym
  set.seed(1)
  v <- c(1, 2, 3) / 10
  A <- matrix(c(0, -v[1], -v[2],
                v[1], 0, -v[3],
                v[2], v[3], 0), 3, 3, byrow=TRUE)
  SSv <- MakeSkewSym(v)

  expect_equal(SSv, A)
  expect_equal(SSv + t(SSv), matrix(0, 3, 3))

  m <- 5
  pts <- seq(0, 1, length.out=m)
  v0 <- function(x) rep(0, length(x))
  v1 <- function(x) x * 2 
  v2 <- function(x) exp(- (x * 2 - 1/4)^2 * 3) - exp(- 3 / 16)
  F1 <- MakeSkewSymFunc(list(v1))
  F2 <- MakeSkewSymFunc(list(v0, v1, v2))
  Fpts1 <- sapply(F1, function(f) f(pts))
  Fpts2 <- sapply(F2, function(f) f(pts))
  Fmat1 <- lapply(seq_len(m), function(i) matrix(Fpts1[i, ], 2, 2))
  Fmat2 <- lapply(seq_len(m), function(i) matrix(Fpts2[i, ], 3, 3))

  # MakeSkewSymFunc
  expect_equal(sapply(Fmat1, `[`, 2, 1), v1(pts))
  expect_equal(sapply(Fmat1, `[`, 1, 2), -v1(pts))
  expect_equal(sapply(Fmat2, `[`, 3, 1), v1(pts))
  expect_equal(sapply(Fmat2, `[`, 1, 3), -v1(pts))
  expect_equal(sapply(Fmat2, `[`, 3, 2), v2(pts))
  expect_equal(sapply(Fmat2, `[`, 2, 3), -v2(pts))
  # Output is skew symmetric
  expect_true(all(sapply(Fmat1, isSkewSym)))
  expect_true(all(sapply(Fmat2, isSkewSym)))

  R1 <- MakeSOMat(v)
  # R2 <- MakeSOMat(c(0, 0, 0))
  R2 <- MakeSOMat(rep(0, 3))
  R3 <- MakeSOMat(rnorm(3))

  # MakeSOMat
  expect_equal(R1, expm::expm(A)) # Make SO(n) matrix works
  expect_equal(R2, diag(3))
  expect_equal(crossprod(MakeSOMat(rnorm(6))), diag(4)) # orthogonal
  expect_equal(crossprod(MakeSOMat(rnorm(10))), diag(5)) # orthonormal

  # DistSO
  expect_equal(DistSO(R1, R3), rotations::rot.dist(rotations::as.SO3(R1), rotations::as.SO3(R3), method='intrinsic') * sqrt(2))
  expect_equal(DistSO(R2, R2), 0)
  # A <- crossprod(R1)
  # expm::logm(A, 'Eigen')
})


test_that('log.SO, exp.SO, distance.SO works', {
  
  mfd <- structure(list(), class='SO')
  x <- c(1, 2, 1) / 5         
  y <- c(-2, 1, -1) / 5
  X <- MakeSOMat(x)
  Y <- MakeSOMat(y)

  # Example from matlab
  z <- matrix(c(1, 1, 0, 0, 0, 2, 0, 0, -1), 3, 3)
  Z <- matrix(c(2.718281828459046, 1.718281828459045, 1.086161269630488, 0, 1, 1.264241117657115, 0, 0, 0.367879441171442), 3, 3)
  A <- cbind(as.numeric(X), as.numeric(Y))

  w <- matrix(c(0, 1, 0, -1, 0, 2, 0, -2, 0) / 5, 3, 3, byrow=TRUE)
  W <- matrix(c(0.980331119030009, 0.193399683419916, 0.039337761939982, -0.193399683419915, 0.901655595150045, 0.386799366839831, 0.039337761939982, -0.386799366839831, 0.921324476120036), 3, 3, byrow=TRUE)

  expect_equal(ExpM(z), Z)
  expect_equal(LogM(Z), z)
  expect_equal(ExpM(w), W)
  expect_equal(ExpMSO3(w), W)
  expect_equal(LogM(W), w)
  expect_equal(LogMSO3(W), w)
  expect_equal(ExpMSO3(matrix(0, 3, 3)), diag(3))
 
  expect_equal(LogM(X), LogMSO3(X))
  expect_equal(LogM(Y), LogMSO3(Y))
  expect_equal(log.SO(mfd, as.numeric(X), as.numeric(X)), 
               matrix(0, 3, 1))
  expect_equal(c(log.SO(mfd, as.numeric(diag(3)), as.numeric(X))), 
               matrix(expm::logm(X), ncol=1)[lower.tri(diag(3))])
  expect_equal(log.SO(mfd, as.numeric(diag(3)), A), unname(cbind(x, y)))

  v1 <- c(-1, -0.5, -1)
  expect_equal(exp.SO(V=as.numeric(v1)), 
               matrix(expm::expm(MakeSkewSym(v1)), ncol=1))
  expect_equal(exp.SO(V=cbind(as.numeric(v1), as.numeric(v1)))[, 1, drop=FALSE], 
               matrix(expm::expm(MakeSkewSym(v1)), ncol=1))

  expect_equal(exp.SO(V=log.SO(X=A)), A)
  expect_equal(exp.SO(mfd, as.numeric(X), log.SO(mfd, as.numeric(X), A)), A)
  expect_equal(log.SO(X=exp.SO(V=log.SO(X=A))), unname(cbind(x, y)))

  expect_equal(rieExp(mfd, V=as.numeric(w[lower.tri(w)])), exp.SO(V=as.numeric(w[lower.tri(w)])))
  expect_equal(rieLog(mfd, X=as.numeric(W)), log.SO(X=as.numeric(W)))

  set.seed(1)
  for (d in 2:3) {
    p0 <- if (d == 2) 1 else if (d == 3) 3
  for (i in 1:10) {
    mu <- exp.SO(V=rnorm(p0))
    Z <- matrix(exp.SO(p=mu, V=rnorm(p0)), d, d)
    expect_true(isSO(Z))
    V <- log.SO(p=mu, X=as.numeric(Z))
    expect_equal(nrow(V), p0)
  }
  }

  # distance.SO
  expect_equal(distance.SO(mfd, as.numeric(X), as.numeric(X)), 0)
  expect_equal(distance.SO(mfd, as.numeric(X), as.numeric(X)), 0)
  expect_equal(distance.SO(mfd, as.numeric(X), as.numeric(diag(3)))^2, sum(x^2))
  expect_equal(distance.SO(mfd, as.numeric(X), as.numeric(Y))^2, 
               sum(log.SO(p=as.numeric(X), X=as.numeric(Y))^2))
})


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
  d <- calcGeomPar.SO(dimTangent = p)
  mfd <- structure(list(), class='SO')
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
  expect_equal(bwMuGCV, 0.1, tolerance=0.1, scale=1)
  expect_equal(bwMuGCVEu, 0.1, tolerance=0.1, scale=1)

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
  expect_equal(max(abs(apply(phi1, 3, function(phi) sum(phi^2)) / m - 1)), 0, scale=1, tolerance=5e-2)
  expect_equal(max(abs(apply(phi2, 3, function(phi) sum(phi^2)) / m - 1)), 0, scale=1, tolerance=5e-2)
  expect_equal(max(abs(apply(phi3, 3, function(phi) sum(phi^2)) / m - 1)), 0, scale=1, tolerance=5e-2)


})


test_that('calcIntDim.SO, calcGeomPar.SO, calcTanDim.SO work', {

  mfd <- structure(1, class='SO')
  for (n in 2:4) {
    X <- diag(n)
    ambient <- length(X)
    intrinsic <- n * (n - 1) / 2
    tangent <- intrinsic

    expect_equal(calcGeomPar(mfd, dimTangent=tangent), n)
    expect_equal(calcIntDim(mfd, dimAmbient=ambient), intrinsic)
    expect_equal(calcIntDim(mfd, dimTangent=tangent), intrinsic)
    expect_equal(calcTanDim(mfd, dimAmbient=ambient), tangent)
    expect_equal(calcTanDim(mfd, dimIntrinsic=intrinsic), tangent)
    expect_error(calcGeomPar(mfd, 1, 2))
    expect_error(calcIntDim(mfd, 1, 2))
    expect_error(calcTanDim(mfd, 1, 2))
  }
})
