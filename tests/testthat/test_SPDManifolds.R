# devtools::load_all()
library(testthat)

test_that('MakePhi.SPD, MakeMfdProcess, GCVFrechetMeanCurve, frechetMeanCurve, RFPCA works for LogEu SPD', {

  mfd <- structure(list(), class=c('LogEu', 'SPD'))
  mfdEu <- structure(1, class='Euclidean')
  n <- 200
  p <- 3
  m <- 11
  pts <- seq(0, 1, length.out=m)
  sparsity <- 2:5
  bwMu <- 0.1 + 1e-5
  kern <- 'epan'
  muScale <- 1
  K <- 5
  lam <- c(0.3, 0.2, 0.1, 0, 0)

  if (p == 1) {
    muList <- list( function(x) x * muScale)
  } else if (p == 3) {
    muList <- list(
      function(x) x * 2 * muScale, 
      function(x) sin(x * 1 * pi) * pi / 2 * 0.6 * muScale,
      function(x) sin(x * 1 * pi) * pi / 2 * 0.6 * muScale,
      function(x) rep(0, length(x))
    )
  } else if (p == 4) {
  }

  set.seed(1)
  d <- calcGeomPar(mfd, dimIntrinsic = p)
  mu <- Makemu(mfd, muList, as.numeric(diag(d)), pts)
  phi <- MakePhi.SPD(mfd, K, mu, pts)
  a <- matrix(phi, m * d^2, K)
  gram <- crossprod(a) / m
  expect_equal(gram, diag(K), tolerance=2 * 1/m)

  samp <- MakeMfdProcess(mfd, n, mu, pts, K = K, lambda=lam, xiFun=rnorm, epsFun=rnorm, sigma2=0.01)
  spSamp <- SparsifyM(samp$XNoisy, samp$T, sparsity)

  yList <- spSamp$Ly
  yListLog <- lapply(yList, ToLower, takeLog=TRUE)
  yListLog0 <- lapply(yList, function(l) {
                       a <- apply(l, 2, function(x) {
                                    m <- logmvec(x, d=d)
                                    diag(m) <- diag(m) / sqrt(2)
                                    m[lower.tri(m, diag=TRUE)]
                            })
  })  

  # ToLower works
  expect_equal(yListLog, yListLog0)

  tList <- spSamp$Lt
  Ymat <- do.call(cbind, yList)
  YmatLog <- do.call(cbind, yListLog)
  Tvec <- do.call(c, tList)
  uniqT <- sort(unique(Tvec))
  tInd <- match(uniqT, pts)
  ord <- order(Tvec)

  # a <- profvis({
  muM <- frechetMeanCurve(mfd, bwMu, kern, xin=Tvec[ord], yin=Ymat[, ord], xout=pts, npoly=1)
  muMEu <- frechetMeanCurve(mfdEu, bwMu, kern, xin=Tvec[ord], yin=YmatLog[, ord], xout=pts, npoly=1)
  muMEuExp <- apply(muMEu, 2, function(x) manifold::ExpM(MakeSym(lowerTri=x, doubleDiag=TRUE)))
  expect_equal(muM, muMEuExp)
  # })
  # bwMuGCV <- GCVFrechetMeanCurve(mfd, yList, tList, kern, 0, 'Sparse')$bOpt
  # bwMuGCVEu <- GCVFrechetMeanCurve(mfdEu, yList, tList, kern, 0, 'Sparse')$bOpt

  res <- RFPCA(yList, tList, list(userBwMu=bwMu, userBwCov=bwMu * 2, kernel=kern, error=TRUE, maxK=K, mfd=mfd))
  resEu <- RFPCA(yListLog, tList, list(userBwMu=bwMu, userBwCov=bwMu * 2, kernel=kern, error=TRUE, maxK=K, mfd=mfdEu, fastEig=FALSE))
  resEu1 <- RFPCA(yListLog, tList, list(userBwMu=bwMu, userBwCov=bwMu * 2, kernel=kern, error=TRUE, maxK=K, mfd=mfdEu, fastEig=TRUE))
  expect_equal(resEu[names(resEu) != 'optns'], resEu1[names(resEu) != 'optns'])

  # The log-Euclidean RFPCA is the same as the RFPCA for the lower triangle of the logged matrices, up to a scale of sqrt(2) 
  expect_equal(res$muObs, apply(resEu$muObs, 2, function(x) expmvec(MakeSym(lowerTri=x, doubleDiag=TRUE), d=2)))
  # plot(res$cov[, , 1, 1], resEu$cov[, , 1, 1]); abline(a=0, b=1/2)
  # plot(res$cov[, , 3, 3], resEu$cov[, , 2, 2]); abline(a=0, b=1)
  expect_equal(unname(pmin(abs(res$xi - resEu$xi * sqrt(2)), abs(res$xi + resEu$xi * sqrt(2)))), matrix(0, n, K))

  a <- matrix(res$phi, ncol=dim(res$phi)[3])
  b <- matrix(aperm(apply(resEu$phi, c(1, 3), function(x) MakeSym(lowerTri=x, doubleDiag=TRUE)), 
                                         c(2, 1, 3)), ncol=dim(res$phi)[3]) / sqrt(2)
  expect_equal(abs(colSums(a * b)), rep(50, K))
  expect_equal(res$sigma2, resEu$sigma2 * 2, tolerance=1e-6)

})



test_that('MakeMfdProcess, GCVFrechetMeanCurve, frechetMeanCurve, RFPCA works for AffInv SPD', {

  mfd <- structure(list(), class=c('AffInv', 'SPD'))
  mfdEu <- structure(1, class='Euclidean')
  set.seed(1)
  n <- 200
  p <- 3
  m <- 11
  pts <- seq(0, 1, length.out=m)
  sparsity <- 2:5
  bwMu <- 0.1 + 1e-5
  kern <- 'epan'
  muScale <- 1
  K <- 5
  lam <- c(0.3, 0.2, 0.1, 0, 0)
  sigma2 <- 0.01

  if (p == 1) {
    muList <- list( function(x) x * muScale)
  } else if (p == 3) {
    muList <- list(
      function(x) x * 2 * muScale, 
      function(x) sin(x * 1 * pi) * pi / 2 * 0.6 * muScale,
      function(x) sin(x * 1 * pi) * pi / 2 * 0.6 * muScale,
      function(x) rep(0, length(x))
    )
  } else if (p == 4) {
  }

  set.seed(1)
  d <- calcGeomPar(mfd, dimIntrinsic = p)
  mu <- Makemu(mfd, muList, as.numeric(diag(d)), pts)
  phi <- MakePhi.SPD(mfd, K, mu, pts)

  samp <- MakeMfdProcess(mfd, n, mu, pts, K = K, lambda=lam, xiFun=rnorm, epsFun=rnorm, sigma2=sigma2)
  spSamp <- SparsifyM(samp$XNoisy, samp$T, sparsity)

  yList <- spSamp$Ly
  yListLog <- lapply(yList, function(l) {
                       a <- apply(l, 2, function(x) {
                                    m <- logmvec(x, d=d)
                                    diag(m) <- diag(m) / sqrt(2)
                                    m[lower.tri(m, diag=TRUE)]
                            })
  })  
  tList <- spSamp$Lt
  Ymat <- do.call(cbind, yList)
  YmatLog <- do.call(cbind, yListLog)
  Tvec <- do.call(c, tList)
  uniqT <- sort(unique(Tvec))
  tInd <- match(uniqT, pts)
  ord <- order(Tvec)

  muM <- frechetMeanCurve(mfd, bwMu, kern, xin=Tvec[ord], yin=Ymat[, ord], xout=pts, npoly=1)
  expect_equal(mu, muM, tolerance=1e-1)

  # LogEu mean
  muMEu <- frechetMeanCurve(mfdEu, bwMu, kern, xin=Tvec[ord], yin=YmatLog[, ord], xout=pts, npoly=1)
  muMEuExp <- apply(muMEu, 2, function(x) manifold::ExpM(MakeSym(lowerTri=x, doubleDiag=TRUE)))
  expect_equal(muM, muMEuExp, tolerance=1e-1)

  res <- RFPCA(yList, tList, 
               list(userBwMu=bwMu, userBwCov=bwMu * 2, 
                    kernel=kern, error=TRUE, maxK=K, mfd=mfd))
  expect_equal(lam, res$lam[seq_len(K)], tolerance=2e-1)
  expect_equal(sigma2, res$sigma2, tolerance=1e-2)

})
