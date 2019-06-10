devtools::load_all()
library(testthat)

test_that('CEScores without error is OK', {
  MC <- 5
  n <- 20
  m <- 30
  K <- 4
  KDirect <- 4
  lambda <- 0.07 ^ (seq_len(K) / 2)
  p <- 3
  basisType <- 'legendre01'
  sparsity <- 5:10
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

  # Generate noiseless samples
  set.seed(1)
  samp <- MakeSphericalProcess(n, mu, pts, K = K, lambda=lambda, basisType=basisType)
  logsamp <- apply(samp$X, 1, function(x) rieLog(structure(list(), class='Sphere'), mu, x))
  samp$X <- aperm(array(logsamp, c(p, m, n)), c(3, 1, 2))
  spSamp <- SparsifyM(samp$X, samp$T, sparsity)

  obsGrid <- sort(unique(unlist(spSamp$Lt)))
  tInd <- match(obsGrid, pts)
  # muObs <- mu[, tInd, drop=FALSE]
  muObs <- matrix(0, p, m)
  trueCov <- samp$phi %*% diag(lambda, length(lambda)) %*% t(samp$phi)
  covObs <- aperm(array(trueCov, c(m, p, m, p)), c(1, 3, 2, 4))[tInd, tInd, , , drop=FALSE]
  phiObs <- array(samp$phi, c(m, p, ncol(samp$phi)))

  res <- CEScores(spSamp$Ly, spSamp$Lt, list(), muObs, obsGrid, covObs, lambda, phiObs, 0)
  est <- t(do.call(cbind, res['xiEst', ]))

  # Estimated scores and true scores are the same
  for (k in seq_len(K)) {
    expect_equal(abs(est[, k] / samp$xi[, k]), rep(1, nrow(est)))
  }

  # # Regularization looks OK
  # res1 <- CEScores(spSamp$Ly, spSamp$Lt, list(), muObs, obsGrid, covObs, lambda, phiObs, 1)
  # est1 <- t(do.call(cbind, res1['xiEst', ]))
  # plot(est1[, 1], samp$xi[, 1])
})
