devtools::load_all()
library(testthat)

n <- 500
p <- 3
e1 <- matrix(1, p, 1)
e2 <- matrix(c(1, rep(0, p - 1)))
m1 <- tcrossprod(e1)
m2 <- tcrossprod(e1, e2)
m3 <- tcrossprod(e2, e1)
m4 <- tcrossprod(e2)
basisMat <- t(sapply(list(m1, m2, m3, m4), as.numeric))

b <- c(1, 1, 2, 1)
trueCov <- m1 + m2 + 2 * m3 + 1 * m4

# Test linComb
expect_equal(linComb(1, matrix(trueCov, nrow=1), makeMat=TRUE), 
             trueCov)
expect_equal(linComb(c(1, 2), 
                     rbind(as.numeric(trueCov), as.numeric(trueCov)), 
                     makeMat=TRUE), 
             3 * trueCov)

# Test smoothRawCov1
sig <- 0.01
xin <- matrix(rnorm(2 * n), n, 2)
beta <- c(1, 2) / 10
yin <- matrix(1, n, 1) %*% matrix(trueCov, nrow=1) + 
  xin %*% matrix(beta, 2, 4)  %*% basisMat + 
  matrix(rnorm(n * 4), n, 4) %*% basisMat
win <- rep(1, n)
expect_equal(smoothRawCov1(xin, yin, win), trueCov, tolerance=0.1)


# # Test smoothCovM
# test_that('smoothCovM agrees with entrywise Lwls2D', {
  # kern <- 'epan'
  # set.seed(1)
  # for (p in 2:3) {
    # # muFun <- function(tvec) {
      # # if (p == 1) {
        # # unname(rbind(tvec))
      # # } else {
        # # unname(rbind(tvec, -tvec, rep(0, p - 2)))
      # # }
    # # }
    # muFun <- function(tvec) {
      # matrix(c(rep(0, p-1), 1), p, length(tvec))
    # }
    # # muFun(1:10)
    # m <- 5:10
    # n <- 30
    # r <- 0.8 # correlation between the first two dimensions
    # lam <- 0.1
    # sig <- 0.05
    # nT <- 20
    # allT <- seq(0, 1, length.out=nT)
    # tList <- lapply(seq_len(n), function(i) {
      # Ni <- sample(c(m, m), 1)
      # sort(sample(allT, Ni))
    # })
    # # Normal around the mean curve. First two dimensions correlated, the rest isotropic. 
    # phi <- function(t) sin(t * 2 * pi) + 1 # Not unit length
    # vList <- lapply(seq_along(tList), function(i) {

                      # # browser()
      # tvec <- tList[[i]]
      # mu <- muFun(tvec)
      # xi <- rnorm(p, sd=sqrt(lam))
      # if (p == 1) {
      # } else if (p > 1) {
        # covMat <- r + diag(1 - r, 2)
        # xi[1:2] <- chol(covMat) %*% xi[1:2]
      # }
      # xi[p] <- 0
      # Xminusmu <- matrix(xi, ncol=1) %*% matrix(phi(tvec), nrow=1)
      # Xminusmu

    # })

    # yList <- lapply(seq_along(tList), function(i) {

      # tvec <- tList[[i]]
      # Xminusmu <- vList[[i]]
      # mfd <- structure(list(), class='Sphere2D')
      # rieExp(mfd, muFun(tvec), Xminusmu)

    # })

    # mu <- muFun(allT)
    # bw <- 0.2
    # resMult <- smoothCovM(yList, tList, mu, allT, bw, kern)
    # # hist(resMult, 20)

    # # Compare with Lwls2D
    # resLwls2D <- smoothCovM2(yList, tList, mu, allT, bw, kern)
    # expect_equal(resMult, resLwls2D)
  # }
# })

# TODO: test that kernel weights sum up to one

test_that('projectCov is correct for spheres', {
  mfd <- structure(1, class='Sphere')
  covR <- array(rnorm(16), c(2, 2, 2, 2))
  covR <- projectCov(mfd, covR, matrix(c(1, 0, 0, 1), 2, 2))
  expect_equal(as.numeric(covR[1, , 1, ]), rep(0, 4))
  expect_equal(as.numeric(covR[2, , 2, ]), rep(0, 4))
  expect_equal(as.numeric(covR[, 1, , 1]), rep(0, 4))
  expect_equal(as.numeric(covR[, 2, , 2]), rep(0, 4))
})


# Some simulation
test_that('smoothCovM and smoothCovM2 estimates the true covariance', {
  MC <- 5
  n <- 50
  m <- 20
  K <- 20
  KDirect <- 4
  lambda <- 0.07 ^ (seq_len(K) / 2)
  p <- 3
  basisType <- 'legendre01'
  sparsity <- 10:15
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
  CreateBasis <- fdapace:::CreateBasis
  samp <- MakeSphericalProcess(n, mu, pts, K = K, lambda=lambda, basisType=basisType)
  spSamp <- SparsifyM(samp$X, samp$T, sparsity)

  trueCov <- samp$phi %*% diag(lambda, length(lambda)) %*% t(samp$phi)

  # smooth
  bw <- 0.2
  kern <- 'epan'
  # res <- smoothCovM(spSamp$Ly, spSamp$Lt, mu, pts, bw, kern)
  res2 <- smoothCovM2(spSamp$Ly, spSamp$Lt, mu, pts, bw, kern)
  projM <- sapply(seq_len(ncol(mu)), 
                  function(i) projectTangent(mfd, mu[, i], projMatOnly=TRUE), 
                  simplify='array')

  for (j1 in seq_len(m)) {
  for (j2 in seq_len(m)) {
    res2[j1, j2, , ] <- projM[, , j1] %*% res2[j1, j2, , ] %*% t(projM[, , j2])
  }
  }

  # resComp <- matrix(aperm(res, c(1, 3, 2, 4)), nrow(trueCov), ncol(trueCov))
  # expect_equal(trueCov, resComp, scale=1, tolerance=5e-2)
  resComp2 <- matrix(aperm(res2, c(1, 3, 2, 4)), nrow(trueCov), ncol(trueCov))
  expect_equal(trueCov, resComp2, scale=1, tolerance=5e-2)
  # e1 <- eigen(trueCov)$vectors[, 1]
  # e1h <- eigen(resComp)$vectors[, 1]
  # matplot(cbind(e1, e1h), type='l')
})
