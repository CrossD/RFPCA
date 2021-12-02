devtools::load_all()
library(testthat)

test_that('MakeInput works', {
  ID <- c('a', 'a', 'a', 'b', 'b')
  tt <- c(1, 2, 3, 1, 2)
  yy <- matrix(rnorm(5 * 2), 5, 2)
  dat <- data.frame(ID, tt, yy, stringsAsFactors=FALSE)
  dat <- dat[sample(nrow(dat)), ]
  res <- MakeInput(dat, 'ID', 'tt')
  res1 <- MakeInput(dat, 'ID', 'tt', c('X1', 'X2'))

  expect_equal(unname(res$Lt), list(c(1, 2, 3), c(1, 2)))
  expect_equal(lapply(unname(res$Ly), unname), 
               list(t(yy[1:3, , drop=FALSE]), 
                    t(yy[4:5, , drop=FALSE])))
  expect_equal(unlist(res$Lid), unique(ID))
  expect_equal(res, res1)
})


test_that('MakeRotMat works', {
  x1 <- c(1, 0, 0)
  x2 <- c(-1, 0, 0)
  x3 <- c(-0.99, sqrt(1 - 0.99^2), 0)

  expect_warning(MakeRotMat(x1, x2), 'Rotation between podal points is arbitrary')

  R2 <- MakeRotMat(x1, x3)

})
