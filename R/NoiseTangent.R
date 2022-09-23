#' Generate noise on the tangent space
#'
#' For each time point (column in mu), generate n isotropic noises distributed as epsFun on the tangent space at mu(t)
#' @param mfd A manifold object created by \code{\link{createM}}
#' @param n Sample size
#' @param mu A matrix or a vector containing the base point(s)
#' @param sigma2 The isotropic variance on the tangent space
#' @param epsFun The distribution for generating the univariate errors with unit variance for the tangent space coordinates
#' 
#' @return An array with size (n, dimIntrinsic, nPts) 
#' @export
NoiseTangent <- function(mfd, n, mu, sigma2=0, epsFun = rnorm) UseMethod('NoiseTangent', mfd)


#' @export
#' @describeIn NoiseTangent Method
NoiseTangent.Sphere <- function(mfd, n, mu, sigma2=0, epsFun = stats::rnorm) {

  m <- ncol(mu)
  dimTangent <- dimAmbient <- nrow(mu)
  res <- vapply(seq_len(m), function(tt) {
           rot <- MakeRotMat(c(rep(0, dimAmbient - 1), 1), mu[, tt])
           eps <- cbind(matrix(epsFun(n * (dimAmbient - 1)), n, dimAmbient - 1), 0) %*% t(rot)
           eps
         }, matrix(0, n, dimTangent)) * sqrt(sigma2)
  # res <- aperm(res, c(1, 3, 2))
  array(res, c(n, dimTangent, m))
}


#' @param sigma2 Noise variance 
#' @export
#' @describeIn NoiseTangent Method
NoiseTangent.Euclidean <- function(mfd, n, mu, sigma2=0, epsFun = stats::rnorm) {

  m <- ncol(mu)
  dimIntrinsic <- dimAmbient <- nrow(mu)
  res <- vapply(seq_len(m), function(tt) {
           eps <- matrix(epsFun(n * dimIntrinsic), n, dimIntrinsic)
           eps
         }, matrix(0, n, dimAmbient)) * sqrt(sigma2)
  # res <- aperm(res, c(1, 3, 2))
  array(res, c(n, dimAmbient, m))
}


#' @export
#' @describeIn NoiseTangent Method
NoiseTangent.SPD <- function(mfd, n, mu, sigma2=0, epsFun = stats::rnorm) {

  m <- ncol(mu)
  dimAmbient <- nrow(mu)
  dimIntrinsic <- calcIntDim(mfd, dimAmbient=dimAmbient)
  res <- array(epsFun(n * dimIntrinsic * m), c(n, dimIntrinsic, m)) * sqrt(sigma2 / 2)
  tmp <- apply(res, c(1, 3), function(x) MakeSym(lowerTri=x, doubleDiag=TRUE))
  res <- aperm(tmp, c(2, 1, 3))

  res

}


NoiseTangent.SO <- function(mfd, n, mu, sigma2=0, epsFun = stats::rnorm) {

  m <- ncol(mu)
  dimAmbient <- nrow(mu)
  dimIntrinsic <- calcIntDim(mfd, dimAmbient=dimAmbient)
  res <- array(epsFun(n * dimIntrinsic * m), c(n, dimIntrinsic, m)) * sqrt(sigma2)
  res

}


# TODO: Either implement this, or replace this with rmfd.*
# NoiseTangent.FlatTorus <- function(mfd, n, mu, sigma2=0, epsFun = stats::rnorm) {

  # m <- ncol(mu)
  # dimIntrinsic <- dimAmbient <- nrow(mu)
  # res <- vapply(seq_len(m), function(tt) {
           # eps <- matrix(epsFun(n * dimIntrinsic), n, dimIntrinsic)
           # eps
         # }, matrix(0, n, dimAmbient)) * sqrt(sigma2)
  # # res <- aperm(res, c(1, 3, 2))
  # array(res, c(n, dimAmbient, m))
# }


# # NoiseTangent.Tree <- function(mfd, n, mu, sigma2=0, epsFun = stats::rnorm) {

  # # m <- ncol(mu)
  # # dimIntrinsic <- dimAmbient <- nrow(mu)
  # # res <- vapply(seq_len(m), function(tt) {
           # # eps <- matrix(epsFun(n * dimIntrinsic), n, dimIntrinsic)
           # # eps
         # # }, matrix(0, n, dimAmbient)) * sqrt(sigma2)
  # # # res <- aperm(res, c(1, 3, 2))
  # # array(res, c(n, dimAmbient, m))
# # }


