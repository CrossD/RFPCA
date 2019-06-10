devtools::load_all()
library(testthat)

test_that('EigenAnalysis works', {
  n <- 500
  p <- 3
  m <- 20
  pts <- seq(0, 1, length.out=m)
  gridSize <- pts[2] - pts[1]
  e1 <- as.numeric(matrix(pts) %*% matrix(1, 1, p))
  e2 <- as.numeric(matrix(pts^2) %*% matrix(c(1, -1, rep(0, p - 2)), nrow=1))

  # Normalize eigenfunctions
  e1 <- e1 / sqrt(sum(e1^2) * gridSize)
  e2 <- e2 / sqrt(sum(e2^2) * gridSize)
  stopifnot(abs(crossprod(e1, e2)) < 1e-15)
  lambda <- c(3, 1)
  trueCov <- aperm(array(tcrossprod(e1) * lambda[1] + tcrossprod(e2) * lambda[2], 
                         dim=c(m, p, m, p)), 
                   c(1, 3, 2, 4))

  res <- EigenAnalysis(trueCov, pts)

  # The dimensions of the input and fitted covariances are the same
  expect_equal(dim(res$fittedCov), dim(trueCov))

  # Eigenvalues are correct
  expect_equal(res$lambda, lambda)

  # Eigenvectors are orthogonal
  expect_equal(crossprod(res$phi) * gridSize, diag(2))

  # Eigenvectors are correct
  expect_equal(abs(crossprod(res$phi, unname(cbind(e1, e2))) * gridSize), diag(2))
})

