# Estimate the isotropic error sigma2
# smooth Whether to smooth the diagonal to obtain the sigma2 estimate (default: TRUE). If FALSE, the mean variance on the diagonal will be used.
EstSigma2 <- function(mfd, yList, tList, mu, covR, regGrid, bwCov, kernel_type, smooth=TRUE) {

  n <- length(yList)
  nT <- ncol(mu) # number of time points
  m <- length(regGrid)
  dimTangent <- nrow(mu)
  dimIntrinsic <- calcIntDim(mfd, dimTangent=dimTangent)

  multyList <- lapply(seq_len(dimTangent), function(j) lapply(yList, `[`, j, ))

  if (smooth) {
    # Smooth the diagonal elements only (t_{ij} = s_{ij})
    diagNoisy <- 
      sapply(seq_len(dimTangent), function(j2) {
      sapply(seq_len(dimTangent), function(j1) {
        rcov <- fdapace:::GetRawCrCovFuncFunc(
          multyList[[j1]], tList, mu[j1, ], 
          multyList[[j2]], tList, mu[j2, ]
          )
        ind <- rcov$tpairn[, 1] == rcov$tpairn[, 2]
        xin <- rcov$tpairn[ind, 1]
        yin <- rcov$rawCCov[ind]
        ord <- order(xin)
        xin <- xin[ord]
        yin <- yin[ord]
        res <- fdapace::Lwls1D(bwCov, kernel_type, xin=xin, yin=yin, xout=regGrid, npoly=0)
        res
      }, simplify='array')
      }, simplify='array')

    varNoisy <- apply(diagNoisy, 1, function(x) sum(diag(as.matrix(x))))

    varEst <- vapply(seq_len(m), 
                     function(tt) sum(diag(matrix(covR[tt, tt, , ], 
                                                  dimTangent, 
                                                  dimTangent))), 
                     0)
    res <- mean(varNoisy - varEst) / dimIntrinsic
  } else {

    # Use the mean of raw covariances to estimate diagonal variance
    totVarNoisy <- 
      sum(sapply(seq_len(dimTangent), function(j) {
        rcov <- fdapace:::GetRawCrCovFuncFunc(
          multyList[[j]], tList, mu[j, ], 
          multyList[[j]], tList, mu[j, ]
          )
        ind <- rcov$tpairn[, 1] == rcov$tpairn[, 2]
        mean(rcov$rawCCov[ind])
      }, simplify='array'))

    totVarEst <- mean(vapply(seq_len(m), 
                      function(tt) sum(diag(matrix(covR[tt, tt, , ], 
                                                   dimTangent, 
                                                   dimTangent))), 
                      0))
    res <- (totVarNoisy - totVarEst) / dimIntrinsic
  }

  if (res < .Machine$double.eps * 100) {
    warning('The error sigma2 is numerically zero or negative. Set to 1e-3 instead')
    res <- 1e-3
  }

  res
}
