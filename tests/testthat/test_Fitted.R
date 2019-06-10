# devtools::load_all()

test_that('PositiveOrth works', {
  M <- t(matrix(c(1, 0, 0,
                  0, -1/sqrt(2), 1/sqrt(2)),
                2, 3, byrow=TRUE))
  Mplus <- t(matrix(c(1, 0, 0,
                      0, 0, 1), 
                    2, 3, byrow=TRUE))
  X <- abind::abind(M, M, along=0)
  expect_equal(PositiveOrth(M), Mplus)
  expect_equal(PositiveOrth(X), unname(abind::abind(Mplus, Mplus, along=0)))
})
