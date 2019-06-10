# devtools::load_all()
# Some simulation
test_that('EstSigma2 works for R^p valued function', {
  MC <- 5
  n <- 500
  m <- 100
  K <- 20
  KDirect <- 4
  lambda <- 0.07 ^ (seq_len(K) / 2)
  p <- 3
  basisType <- 'legendre01'
  sparsity <- 20
  sigma2 <- 1e-1
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

  # Generate noisy samples
  set.seed(1)
  samp <- MakeSphericalProcess(n, mu, pts, K = K, lambda=lambda, basisType=basisType)
  samp$X <- samp$X + rnorm(n * m * p, sd=sqrt(sigma2))
  spSamp <- SparsifyM(samp$X, samp$T, sparsity)
  # matplot(t(samp$X[1, , ]), type='l')

  trueCov <- array(samp$phi %*% diag(lambda, length(lambda)) %*% t(samp$phi),
                   c(m, p, m, p))
  trueCov <- aperm(trueCov, c(1, 3, 2, 4))

  # smooth
  bw <- 0.1
  kern <- 'epan'
  sigma2R <- EstSigma2(mfd, spSamp$Ly, spSamp$Lt, mu, trueCov, pts, bw, kern)
  expect_equal(sigma2R, sigma2, tolerance=1e-1, scale=1)
})
