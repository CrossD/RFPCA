MakePhi.SPD <- 
  function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           type=c('cos', 'sin', 'fourier', 'legendre01', 'poly')) {

  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }
  dimAmbient <- nrow(mu)
  dimIntrinsic <- calcIntDim(mfd, dimAmbient = dimAmbient)
  dimTangent <- calcTanDim(mfd, dimAmbient = dimAmbient)
  dimMat <- round(sqrt(dimAmbient))
  if (dimAmbient <= 1) {
    stop('mu must have more than 1 rows')
  }
  m <- length(pts)
  if (ncol(mu) != m) {
    stop('mu must have the same number of columns as the length of pts')
  }

  if (min(pts) < 0 || max(pts) > 1) {
    stop('The range of pts should be within [0, 1]')
  }

  ptsSqueeze <- do.call(c, lapply(seq_len(dimIntrinsic), function(i) {
                          pts + (max(pts) + mean(diff(pts))) * (i - 1)
                        }))
  ptsSqueeze <- ptsSqueeze / max(ptsSqueeze)
  
  phi <- array(fdapace::CreateBasis(K, ptsSqueeze, type), c(m, dimIntrinsic, K)) / 
    sqrt(dimIntrinsic)
  # Generate the lower triangle, and then symmetrize
  tmp <- apply(phi, c(1, 3), function(x) {
    x <- MakeSym(lowerTri = x, doubleDiag=TRUE)
    x 
  }) / sqrt(2)
  phi <- aperm(tmp, c(2, 1, 3))

  dimnames(phi) <- list(t=pts, j=seq_len(dimTangent), k=seq_len(K))
  phi
}


# Use a trick to convert a matrix representation to the log-lower trianglar representation
# l A matrix whose columns are vectorized matrices
ToLower <- function(l, takeLog=TRUE) {
  l <- as.matrix(l)
  d <- round(sqrt(nrow(l)))
  a <- apply(l, 2, function(x) {
               if (takeLog) {
                 m <- logmvec(x, d=d)
               } else {
                 m <- matrix(x, d, d)
               }
               diag(m) <- diag(m) / sqrt(2)
               m[lower.tri(m, diag=TRUE)]
  })

  a
}

# from {manifold}
logmvec <- function(x, d) manifold::LogM(matrix(x, d, d))
expmvec <- function(x, d) manifold::ExpM(matrix(x, d, d))


