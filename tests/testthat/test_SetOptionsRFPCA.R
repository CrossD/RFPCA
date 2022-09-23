# devtools::load_all()
library(testthat)

test_that('AppendOptionsRFPCA  works', {

  m <- 1
  p <- 40
  tpts <- seq(0, 1, length.out=m)
  spts <- seq(0, 1, length.out=p)
  mfdHS <- structure(1, class=c('HS', 'L2')) # Square-root densities
  K <- 15
  lambda <- 0.05 ^ (seq_len(K) / 3)
  a <- 1
  b <- 1
  mu <- matrix(dbeta(spts, b, a) + 1)
  mu <- sqrt(mu / mean(mu) * (p - 1) / p)
  set.seed(1)
  samp <- MakeMfdProcess(mfdHS, 10, mu, tpts, K = K, lambda=lambda)
  spSamp <- SparsifyM(samp$X, samp$T, m)
  yList <- spSamp$Ly
  tList <- spSamp$Lt
  opt <- SetOptionsRFPCA(yList, tList, optns=list())

  newOpt <- list(maxK = K, FVEthreshold=0.95)
  res <- AppendOptionsRFPCA(opt, maxK = K, FVEthreshold=0.95)

  expect_equal(res$maxK, K)
  expect_equal(res$FVEthreshold, 0.95)
})
