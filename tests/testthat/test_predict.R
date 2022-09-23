# devtools::load_all()
library(testthat)

# Some simulation
test_that('predict.RFPCA works', {
  nTest <- 20
  nTrain <- 200
  n <- nTrain + nTest
  indTrain <- seq_len(nTrain)
  indTest <- seq(nTrain + 1, nTrain + nTest)
  m <- 20
  K <- 5
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

  samp <- MakeSphericalProcess(n, mu, pts, K = K, lambda=lambda, basisType=basisType, sigma2=sigma2)
  spSamp <- SparsifyM(samp$X, samp$T, sparsity)
  yList <- spSamp$Ly[indTrain]
  tList <- spSamp$Lt[indTrain]
  yListTest <- spSamp$Ly[indTest]
  tListTest <- spSamp$Lt[indTest]

  # smooth
  bw <- 0.2
  kern <- 'epan'
  assumeError <- TRUE

  set.seed(1)
  resSp <- RFPCA(yList, tList, list(userBwMu=bw, userBwCov=bw * 2, kernel=kern, error=assumeError, maxK=K))


  KTrain <- K
  # In sample
  insamp <- predict(resSp, yList, tList, type='xi')
  expect_equal(abs(sapply(seq_len(KTrain), 
                          function(k) cor(insamp[, k], resSp$xi[indTrain, k]))), 
               rep(1, KTrain))

  KTest <- 3
  # Out of sample
  pred <- predict(resSp, yListTest, tListTest, type='xi')
  expect_equal(abs(sapply(seq_len(KTest), function(k) cor(pred[, k], samp$xi[indTest, k]))), rep(1, KTest), tolerance=0.1)
})
