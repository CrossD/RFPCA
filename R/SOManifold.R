MakePhi.SO <- 
  function(mfd, K, mu, pts=seq(0, 1, length.out=ncol(mu)), 
           type=c('cos', 'sin', 'fourier', 'legendre01', 'poly')) {

  type <- match.arg(type)
  if (!is.matrix(mu)) {
    stop('mu has to be a matrix')
  }
  dimAmbient <- nrow(mu)
  dimIntrinsic <- calcIntDim.SO(mfd, dimAmbient=dimAmbient)
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

  dimnames(phi) <- list(t=pts, j=seq_len(dimIntrinsic), k=seq_len(K))
  phi

  # vapply(seq_len(m), function(tt) {
    # vapply(seq_len(K), function(k) {
      # browser()
      # MakeSkewSym(epsFun)
    # }, rep(0, dimAmbient))
  # }, matrix(dimAmbient, K))
  # dimnames(phi1) <- list(t=pts, j=seq_len(dimAmbient), k=seq_len(K))

  # # Rotate to the locations of mu
  # basePt <- c(rep(0, dimTangent), 1)
  # phi <- sapply(seq_along(pts), function(i) {
    # R <- MakeRotMat(basePt, mu[, i])
    # R %*% phi1[i, , , drop=TRUE]
  # }, simplify=FALSE)
  # phi <- do.call(abind::abind, list(phi, along=0))
  # dimnames(phi) <- dimnames(phi1)

  # phi
}
