# devtools::load_all()
# library(testthat)

# test_that('Projection and projection matrix works', {

  # n <- 200
  # p <- 4
  # m <- 50
  # pts <- seq(0, 1, length.out=m)
  # sparsity <- 2:5
  # bwMu <- 0.1
  # kern <- 'epan'

  # if (p == 3) {
    # muList <- list(
      # function(x) x * 2 * muScale, 
      # function(x) sin(x * 1 * pi) * pi / 2 * 0.6 * muScale,
      # function(x) rep(0, length(x))
    # )
  # } else if (p == 4) {
    # muList <- list(
      # function(x) x * sqrt(2), 
      # function(x) x * sqrt(2), 
      # function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
      # function(x) rep(0, length(x))
    # )
  # }

  # mfd <- structure(list(), class='Sphere')
  # mu <- Makemu(mfd, muList, c(rep(0, p - 1), 1), pts)
  # samp <- MakeSphericalProcess(n, mu, pts, K = 3, lambda=c(0.3, 0.2, 0.1), xiFun=rnorm, epsFun=rnorm, sigma2=0.01)
  # spSamp <- SparsifyM(samp$XNoisy, samp$T, sparsity)

  # yList <- spSamp$Ly
  # tList <- spSamp$Lt
  # Ymat <- do.call(cbind, yList)
  # Tvec <- do.call(c, tList)
  # uniqT <- sort(unique(Tvec))
  # tInd <- match(uniqT, pts)
  # ord <- order(Tvec)

  # # a <- profvis({
  # muM <- frechetMeanCurve(mfd, bwMu, kern, xin=Tvec[ord], yin=Ymat[, ord], xout=pts, npoly=0)
  # # })
  # bwMuGCV <- GCVFrechetMeanCurve(mfd, yList, tList, kern, 0, 'Sparse')$bOpt
  # bwMuGCVEu <- GCVFrechetMeanCurve(structure(1, class='Euclidean'), yList, tList, kern, 0, 'Sparse')$bOpt
  # matplot(t(mu), type='l', lty=1)
  # matplot(t(muM), type='l', lty=2)
  # expect_equal(mu, muM, tolerance=2e-1)
  # expect_equal(bwMuGCV, 0.1, tolerance=0.05, scale=1)
  # expect_equal(bwMuGCVEu, 0.1, tolerance=0.05, scale=1)

# })
