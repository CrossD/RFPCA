# Bootstrap two-sample test
# samp A list of lists containing the Riemannian functional data
# group A vector of indices for the two populations
# B The number of bootstrap samples
# mu0 The base point of the logarithm map used for tau
BootTest2 <- function(samp, group, mu0, 
                      mfdCmp=structure(1, class=c('HS', 'L2')), 
                      method=c('L2', 'proj', 'all'),
                      K, FVE, 
                      B=500, paired=FALSE, 
                      parallel = c("no", "multicore", "snow"), 
                      ncpus=1L, 
                      RFPCAoptns=list()) {

  method <- match.arg(method)
  parallel <- match.arg(parallel)
  mfd <- RFPCAoptns[['mfd']]
  group <- factor(group)
  uniGroups <- sort(unique(group))
  stopifnot(length(uniGroups) == 2)
  allInd1 <- which(group == uniGroups[1])
  allInd2 <- which(group == uniGroups[2])
  n1 <- length(allInd1)
  n2 <- length(allInd2)
  n <- length(group)
  ns <- nrow(mu0)

  # resBoth <- RFPCA(samp$Ly, samp$Lt, RFPCAoptns)
  samp1All <- list(Ly = samp$Ly[allInd1], Lt=samp$Lt[allInd1])
  samp2All <- list(Ly = samp$Ly[allInd2], Lt=samp$Lt[allInd2])
  res1 <- RFPCA(samp1All$Ly, samp1All$Lt, RFPCAoptns)
  res2 <- RFPCA(samp2All$Ly, samp2All$Lt, RFPCAoptns)
  mu1Hat <- project(mfdCmp, res1$muObs)
  mu2Hat <- project(mfdCmp, res2$muObs)

      # Determine K
  if (method %in% c('all', 'proj')) {
    if (missing(K)) {
      if (missing(FVE)) {
        stop('K and FVE cannot both be missing')
      } else {
        asympCov1 <- GetAsympCov(samp1All, mu1Hat, mu0, mfdCmp)
        asympCov2 <- GetAsympCov(samp2All, mu2Hat, mu0, mfdCmp)
        asympCov <- n / n1 * asympCov1 + n / n2 * asympCov2
        eig <- eigen(asympCov)
        ev <- eig$values / (ns - 1)
        K <- vapply(FVE, function(fve) {
          min(which(cumsum(ev) / sum(ev) >= fve))
        }, 1)
      }
    }
  }

  # Use mu0 to form chart tau. The default is the logarithm map based at mu0
  # mu0 <- project(mfdCmp, resBoth$muObs[, 1])
  tau <- function(X) rieLog(mfdCmp, mu0, X)
  mfdL2 <- structure(1, class='L2')

  # t1, t2: observed statistics
  # t0_1, t0_2: appropriate centers of the observed statistics
  TL2 <- function(t1, t2, t0_1, t0_2) {
    n * norm(mfdL2, U=((t1 - t0_1) - (t2 - t0_2)))^2 
  }

  TL2Dist <- function(t1, t2, t0_1, t0_2) {
    rmat <- MakeRotMat(Normalize(t0_2), Normalize(t0_1))
    distance(mfdCmp, t1, rmat %*% t2)
  }

  TProj <- function(t1, t2, t0_1, t0_2, phi, lam) {
    projs <- metric(mfdL2, 
                    U=(t1 - t0_1) - (t2 - t0_2), 
                    V=phi)
    unname(n * sum(projs^2 / lam))
  }

  if (paired) {
    dat <- seq_len(n / 2)
  } else if (!paired) {
    dat <- seq_len(n)
  }

  if (method == 'L2') {
    # L2 method

    if (!paired) {
      statistic <- function(dat, ind, ...) {
        ddd <- list(...)
        muCmp1 <- ddd$muCmp1
        muCmp2 <- ddd$muCmp2
        ind1 <- ind[ind %in% allInd1]
        ind2 <- ind[ind %in% allInd2]
        res1 <- RFPCA(samp$Ly[ind1], samp$Lt[ind1], RFPCAoptns)
        res2 <- RFPCA(samp$Ly[ind2], samp$Lt[ind2], RFPCAoptns)
        mu1star <- project(mfdCmp, res1$muObs)
        mu2star <- project(mfdCmp, res2$muObs)
        T <- c(TL2(tau(mu1star), tau(mu2star), 
                   tau(muCmp1), tau(muCmp2)),
               TL2Dist(mu1star, mu2star, muCmp1, muCmp2))
        T <- stats::setNames(T, c('L2', 'Dist'))
        T
      }
      # debug(statistic)

    } else if (paired) {
      statistic <- function(dat, ind, ...) {
        ddd <- list(...)
        muCmp1 <- ddd$muCmp1
        muCmp2 <- ddd$muCmp2
        ind1 <- allInd1[ind]
        ind2 <- allInd2[ind]
        res1 <- RFPCA(samp$Ly[ind1], samp$Lt[ind1], RFPCAoptns)
        res2 <- RFPCA(samp$Ly[ind2], samp$Lt[ind2], RFPCAoptns)
        mu1star <- project(mfdCmp, res1$muObs)
        mu2star <- project(mfdCmp, res2$muObs)
        T <- c(TL2(tau(mu1star), tau(mu2star), 
                   tau(muCmp1), tau(muCmp2)),
               TL2Dist(mu1star, mu2star, muCmp1, muCmp2))
        T <- stats::setNames(T, c('L2', 'Dist'))
        T
      }
    }
  } else if (method == 'proj') {
    # Projection method

    if (!paired) {

      statistic <- function(dat, ind, ...) {
        ddd <- list(...)
        muCmp1 <- ddd$muCmp1
        muCmp2 <- ddd$muCmp2

        ind1 <- ind[ind %in% allInd1]
        ind2 <- ind[ind %in% allInd2]
        samp1 <- list(Ly = samp$Ly[ind1], Lt = samp$Lt[ind1])
        samp2 <- list(Ly = samp$Ly[ind2], Lt = samp$Lt[ind2])

        res1 <- RFPCA(samp1$Ly, samp1$Lt, RFPCAoptns)
        res2 <- RFPCA(samp2$Ly, samp2$Lt, RFPCAoptns)
        mu1star <- project(mfdCmp, res1$muObs)
        mu2star <- project(mfdCmp, res2$muObs)

        asympCov1 <- GetAsympCov(samp1, mu1star, mu0, mfdCmp)
        asympCov2 <- GetAsympCov(samp2, mu2star, mu0, mfdCmp)
        asympCov <- n / n1 * asympCov1 + n / n2 * asympCov2
        eig <- eigen(asympCov)
        ev <- eig$values / (ns - 1)
        ef <- eig$vectors * sqrt(ns - 1)

        T <- vapply(K, function(k) {
          TProj(tau(mu1star), tau(mu2star), 
                tau(muCmp1), tau(muCmp2), 
                ef[, seq_len(k), drop=FALSE], 
                ev[seq_len(k)])
        }, 1)
        stats::setNames(T, K)
      }

    } else if (paired) {
      stop()
    }

  } else if (method == 'all') {

    if (!paired) {

      statistic <- function(dat, ind, ...) {
        ddd <- list(...)
        muCmp1 <- ddd$muCmp1
        muCmp2 <- ddd$muCmp2
        # tau <- ddd$tauCmp

        ind1 <- ind[ind %in% allInd1]
        ind2 <- ind[ind %in% allInd2]
        sampNew1 <- list(Ly = samp$Ly[ind1], Lt = samp$Lt[ind1])
        sampNew2 <- list(Ly = samp$Ly[ind2], Lt = samp$Lt[ind2])

        res1 <- RFPCA(sampNew1$Ly, sampNew1$Lt, RFPCAoptns)
        res2 <- RFPCA(sampNew2$Ly, sampNew2$Lt, RFPCAoptns)
        mu1star <- project(mfdCmp, res1$muObs)
        mu2star <- project(mfdCmp, res2$muObs)

        asympCov1 <- GetAsympCov(sampNew1, mu1star, mu0, mfdCmp)
        asympCov2 <- GetAsympCov(sampNew2, mu2star, mu0, mfdCmp)
        asympCov <- n / n1 * asympCov1 + n / n2 * asympCov2
        eig <- eigen(asympCov)
        ev <- eig$values / (ns - 1)
        ef <- eig$vectors * sqrt(ns - 1)

        tL2 <- TL2(tau(mu1star), tau(mu2star), 
                    tau(muCmp1), tau(muCmp2))
        tDist <- TL2Dist(mu1star, mu2star, muCmp1, muCmp2)
        tProj <- vapply(K, function(k) {
          TProj(tau(mu1star), tau(mu2star), 
                tau(muCmp1), tau(muCmp2), 
            ef[, seq_len(k), drop=FALSE], 
            ev[seq_len(k)])
        }, 1)
        T <- stats::setNames(c(tL2, tDist, tProj), c('L2', 'Dist', K))
        T
      }

    } else if (paired) {
      stop()
    }

  }

  if (!paired) {
    res <- boot::boot(dat, statistic, B, strata=group, parallel=parallel, ncpus=ncpus, 
                      muCmp1=project(mfdCmp, res1$muObs), muCmp2=project(mfdCmp, res2$muObs))
    res$t0 <- statistic(dat, seq_len(n), muCmp1=mu0, muCmp2=mu0)
  } else if (paired) {
    res <- boot::boot(dat, statistic, B, parallel=parallel, ncpus=ncpus,
                      muCmp1=project(mfdCmp, res1$muObs), muCmp2=project(mfdCmp, res2$muObs))
    res$t0 <- statistic(dat, seq_len(n / 2), muCmp1=mu0, muCmp2=mu0)
  }

  rbind(T = res$t0, 
        pv = unname(colMeans(res$t > 
                             matrix(res$t0, B, length(res$t0), byrow = TRUE))))
}


# BootTest2_dist <- function(samp, group, mfdCmp=structure(1, class=c('HS', 'L2')), B=500, paired=FALSE, parallel = c("no", "multicore", "snow"), ncpus=1L, RFPCAoptns=list()) {

  # parallel <- match.arg(parallel)
  # mfd <- RFPCAoptns[['mfd']]
  # group <- factor(group)
  # uniGroups <- sort(unique(group))
  # stopifnot(length(uniGroups) == 2)
  # allInd1 <- which(group == uniGroups[1])
  # allInd2 <- which(group == uniGroups[2])
  # n1 <- length(allInd1)
  # n2 <- length(allInd2)
  # n <- length(group)
  # # resBoth <- RFPCA(samp$Ly, samp$Lt, RFPCAoptns)
  # res1 <- RFPCA(samp$Ly[allInd1], samp$Lt[allInd1], RFPCAoptns)
  # res2 <- RFPCA(samp$Ly[allInd2], samp$Lt[allInd2], RFPCAoptns)
  # mu1Hat <- project(mfdCmp, res1$muObs)
  # mu2Hat <- project(mfdCmp, res2$muObs)

  # rmat <- MakeRotMat(Normalize(mu2Hat), Normalize(mu1Hat))

  # T0 <- distance(mfdCmp, mu1Hat, mu2Hat)
  # Tstar <- function(mu1, mu2) {
    # distance(mfdCmp, mu1, rmat %*% mu2)
  # }
  # # Tstar(tau(res1$muObs), tau(res2$muObs))

  # if (!paired) {
    # dat <- seq_len(n)
    # statistic <- function(dat, ind) {
      # ind1 <- ind[ind %in% allInd1]
      # ind2 <- ind[ind %in% allInd2]
      # res1 <- RFPCA(samp$Ly[ind1], samp$Lt[ind1], RFPCAoptns)
      # res2 <- RFPCA(samp$Ly[ind2], samp$Lt[ind2], RFPCAoptns)
      # mu1star <- project(mfdCmp, res1$muObs)
      # mu2star <- project(mfdCmp, res2$muObs)
      # Tstar(mu1star, mu2star)
    # }
    # # debug(statistic)

    # res <- boot::boot(dat, statistic, B, strata=group, parallel=parallel, ncpus=ncpus)
  # } else if (paired) {
    # dat <- seq_len(n / 2)
    # statistic <- function(dat, ind) {
      # ind1 <- allInd1[ind]
      # ind2 <- allInd2[ind]
      # res1 <- RFPCA(samp$Ly[ind1], samp$Lt[ind1], RFPCAoptns)
      # res2 <- RFPCA(samp$Ly[ind2], samp$Lt[ind2], RFPCAoptns)
      # mu1star <- project(mfdCmp, res1$muObs)
      # mu2star <- project(mfdCmp, res2$muObs)
      # Tstar(mu1star, mu2star)
    # }
    # res <- boot::boot(dat, statistic, B, parallel=parallel, ncpus=ncpus)
 # }

  # res$t0 <- T0
  # res
# }



# Bootstrap one-sample test
# samp A list of lists containing the Riemannian functional data
# B The number of bootstrap samples
BootTest1 <- function(samp, mu0, mfdCmp=structure(1, class=c('HS', 'L2')), 
                      method=c('L2', 'proj', 'all'),
                      K, FVE, 
                      B=500, parallel = c("no", "multicore", "snow"), 
                      ncpus=1L, RFPCAoptns=list()) {

  method <- match.arg(method)
  mu0 <- matrix(mu0)
  parallel <- match.arg(parallel)
  if (ncpus > 1L && parallel == 'no') {
    stop('Specify `parallel` since you set `ncpus`')
  }
  n <- length(samp$Ly)
  ns <- nrow(mu0)
  # if (method == 'proj') {
    # if (missing(K) && missing(FVE)) {
      # stop('K and FVE cannot be both missing for the projection method')
    # }
    # if (!missing(K)) {
      # RFPCAoptns <- AppendOptionsRFPCA(RFPCAoptns, maxK=K)
    # } else if (!missing(FVE)) {
      # RFPCAoptns <- AppendOptionsRFPCA(RFPCAoptns, FVEthreshold=FVE)
    # }
  # }
  res <- RFPCA(samp$Ly, samp$Lt, RFPCAoptns)
  stopifnot(all(dim(res$muObs) == dim(mu0)))

  if (method %in% c('proj', 'all')) {
    if (missing(K)) {
      if (missing(FVE)) {
        stop('K and FVE cannot both be missing')
      } else {
        asympCov <- GetAsympCov(samp, project(mfdCmp, res$muObs), mu0, mfdCmp)
        eig <- eigen(asympCov)
        ev <- eig$values / (ns - 1)
        K <- vapply(FVE, function(fve) {
                      min(which(cumsum(ev) / sum(ev) >= fve))
                      }, 1)
      }
    }
    mfdL2 <- structure(1, class='L2')
  }

  dat <- seq_len(n)
  tau <- function(X) rieLog(mfdCmp, mu0, X)
  tauCmp <- function(X) rieLog(mfdCmp, project(mfdCmp, res$muObs), X)

  if (method == 'L2') {
    T <- function(mu1, mu2) {
      mu1 <- project(mfdCmp, mu1)
      mu2 <- project(mfdCmp, mu2)
      as.numeric(distance(mfdCmp, mu1, mu2))
    }

    statistic <- function(dat, ind, ...) {
      ddd <- list(...)
      muCmp <- ddd$muCmp
      resBoot <- RFPCA(samp$Ly[ind], samp$Lt[ind], RFPCAoptns)
      T(resBoot$muObs, muCmp)
    }
  } else if (method == 'proj') {

    statistic <- function(dat, ind, ...) {
      ddd <- list(...)
      muCmp <- ddd$muCmp
      tau <- ddd$tauCmp
      sampNew <- list(Ly=samp$Ly[ind], Lt=samp$Lt[ind])
      resBoot <- RFPCA(sampNew$Ly, sampNew$Lt, RFPCAoptns)
      asympCov <- GetAsympCov(sampNew, project(mfdCmp, resBoot$muObs), muCmp, mfdCmp)
      eig <- eigen(asympCov)
      ev <- eig$values / (ns - 1)
      ef <- eig$vectors * sqrt(ns - 1)
      T <- vapply(K, function(k) {
        projs <- sqrt(n) * 
          metric(mfdL2, 
                 U=tau(project(mfdCmp, resBoot$muObs)) - tau(project(mfdCmp, muCmp)), 
                 V=ef[, seq_len(k), drop=FALSE])
        unname(sum(projs^2 / ev[seq_len(k)]))
      }, 1)
      T
    }
 
  } else if (method == 'all') {

    Tdist <- function(mu1, mu2) {
      mu1 <- project(mfdCmp, mu1)
      mu2 <- project(mfdCmp, mu2)
      as.numeric(distance(mfdCmp, mu1, mu2))
    }

    statistic <- function(dat, ind, ...) {
      ddd <- list(...)
      muCmp <- project(mfdCmp, ddd$muCmp)
      tau <- ddd$tauCmp
      sampNew <- list(Ly=samp$Ly[ind], Lt=samp$Lt[ind])
      resBoot <- RFPCA(sampNew$Ly, sampNew$Lt, RFPCAoptns)
      asympCov <- GetAsympCov(sampNew, project(mfdCmp, resBoot$muObs), muCmp, mfdCmp)
      eig <- eigen(asympCov)
      ev <- eig$values / (ns - 1)
      ef <- eig$vectors * sqrt(ns - 1)
      Tproj <- vapply(K, function(k) {
        projs <- sqrt(n) * 
          metric(mfdL2, 
                 U=tau(project(mfdCmp, resBoot$muObs)) - tau(project(mfdCmp, muCmp)), 
                 V=ef[, seq_len(k), drop=FALSE])
        unname(sum(projs^2 / ev[seq_len(k)]))
      }, 1)
      T <- stats::setNames(c(Tdist(resBoot$muObs, muCmp), Tproj),
                    c('dist', K))
      T
    }
  }

  T0 <- statistic(dat, seq_len(n), muCmp=mu0, tauCmp=tau)
  res <- boot::boot(dat, statistic, B, parallel=parallel, ncpus=ncpus, 
                    muCmp=project(mfdCmp, res$muObs), tauCmp=tauCmp)
  res$t0 <- T0

  rbind(T = T0, 
        pv = unname(colMeans(res$t > 
                             matrix(res$t0, B, length(res$t0), byrow = TRUE))))
}


AsympTest2 <- function(samp, group, mu0, mfdCmp=structure(1, class=c('HS', 'L2')), 
                       paired=FALSE, RFPCAoptns=list(), 
                       method=c('L2', 'proj'),
                       K, FVE, 
                       pvMethod=c('MC', 'bx'), pvMC=1e5) {

  method <- match.arg(method)
  pvMethod <- match.arg(pvMethod)
  mu0 <- matrix(mu0)
  group <- factor(group)
  uniGroups <- sort(unique(group))
  stopifnot(length(uniGroups) == 2)
  allInd1 <- which(group == uniGroups[1])
  allInd2 <- which(group == uniGroups[2])
  n1 <- length(allInd1)
  n2 <- length(allInd2)
  n <- n1 + n2
  ns <- length(mu0)
  samp1 <- list(Ly = samp$Ly[allInd1], Lt = samp$Lt[allInd1])
  samp2 <- list(Ly = samp$Ly[allInd2], Lt = samp$Lt[allInd2])

  # resBoth <- RFPCA(samp$Ly, samp$Lt, RFPCAoptns)
  res1 <- RFPCA(samp1$Ly, samp1$Lt, RFPCAoptns)
  res2 <- RFPCA(samp2$Ly, samp2$Lt, RFPCAoptns)
  # mu1Hat <- project(mfdCmp, res1$muObs)
  # mu2Hat <- project(mfdCmp, res2$muObs)
  mu1Hat <- res1$muObs
  mu2Hat <- res2$muObs

  tau <- function(X) rieLog(mfdCmp, mu0, X)

  # # Distance based 
  # rmat <- MakeRotMat(Normalize(mu2Hat), Normalize(mu1Hat))
  # T0_dist <- distance(mfdCmp, mu1Hat, mu2Hat)
  # Tstar_dist <- function(mu1, mu2) {
    # distance(mfdCmp, mu1, rmat %*% mu2)
  # }

    if (!paired) {
      # Get the distribution of the L2 norm of a Gaussian rv with cov equal to asympCov
      asympCov1 <- GetAsympCov(samp1, mu1Hat, mu0, mfdCmp)
      asympCov2 <- GetAsympCov(samp2, mu2Hat, mu0, mfdCmp)
      asympCov <- n / n1 * asympCov1 + n / n2 * asympCov2
      eig <- eigen(asympCov)
      ev <- eig$values / (ns - 1)
      ef <- eig$vectors * sqrt(ns - 1)

      # ev1 <- eigen(asympCov1)$values / (ns - 1)
      # ev2 <- eigen(asympCov2)$values / (ns - 1)
      # ev <- c(ev1[ev1 > 1e-8] * n / n1, ev2[ev2 > 1e-8] * n / n2)

      if (method == 'L2') {
        T <- unname(n * norm(mfdCmp, U=tau(mu1Hat) - tau(mu2Hat))^2)

        ev <- ev[ev > 1e-8]
        if (pvMethod == 'MC') {
          pv <- pvSim(T, ev, pvMC)
        } else {
          pv <- dr::dr.pvalue(ev[ev > 0], T, chi2approx='bx')[['pval.adj']] # Not very accurate
        }
        res <-   c(T=T, pv=pv)

      } else if (method == 'proj') {

        if (missing(K)) {
          if (missing(FVE)) {
            stop('K and FVE cannot both be missing')
          } else {
            K <- vapply(FVE, function(fve) {
              min(which(cumsum(ev) / sum(ev) >= fve))
            }, 1)
          }
        }
        mfdL2 <- structure(1, class='L2')
        res <- vapply(K, function(k) {
          projs <- sqrt(n) * 
            metric(mfdL2, 
                   U=tau(mu1Hat) - tau(mu2Hat), 
                   V=ef[, seq_len(k), drop=FALSE])
          T <- unname(sum(projs^2 / ev[seq_len(k)]))
          pv <- stats::pchisq(T, k, lower.tail=FALSE)
          c(T=T, pv=pv)
        }, c(1, 1))
      }
  } else if (paired) {
    stop()
  }
  
  res
}


# #' method The method for the asymptotic test statistic. Either 'L2', the scaled L2 distance, or 'proj', the scaled sum of the first K projections onto the eigenfunctions of the asymptotic covariance of tau(hat(mu))
# #' K The number of components to project onto. 
# #' FVE An alternative way to specify the number of components to project onto. `K` and `FVE` cannot be both missing if method == 'proj'
# #' pvMethod Ignored if method == 'proj'. Otherwise, obtain the p-value of a weighted sum of chi-square_1 r.v.s by Monte carlo ('MC') or Bentler-Xie method
AsympTest1 <- function(samp, mu0, mfdCmp=structure(1, class=c('HS', 'L2')),
                       RFPCAoptns=list(), 
                       method=c('L2', 'proj'),
                       K, FVE, 
                       pvMethod=c('MC', 'bx'), 
                       pvMC=1e5) {

  method <- match.arg(method)
  pvMethod <- match.arg(pvMethod)
  mu0 <- matrix(mu0)
  ns <- length(mu0)
  n <- length(samp$Ly)
  # if (method == 'proj') {
    # if (missing(K) && missing(FVE)) {
      # stop('K and FVE cannot be both missing for the projection method')
    # }
    # if (!missing(K)) {
      # RFPCAoptns <- AppendOptionsRFPCA(RFPCAoptns, maxK=K)
    # } else if (!missing(FVE)) {
      # RFPCAoptns <- AppendOptionsRFPCA(RFPCAoptns, FVEthreshold=FVE)
    # }
  # }
  res <- RFPCA(samp$Ly, samp$Lt, RFPCAoptns)
  muHat <- res[['muObs']]

  asympCov <- GetAsympCov(samp, muHat, mu0, mfdCmp)
  # Get the distribution of the L2 norm of a Gaussian rv with cov equal to asympCov
  # The squared norm of a Gaussian process is a sum of chisq variables weighted by the eigenvalues
  eig <- eigen(asympCov)
  ev <- eig$values / (ns - 1)
  ef <- eig$vectors * sqrt(ns - 1)

  tau <- function(X) rieLog(mfdCmp, mu0, X)

  if (method == 'L2') {
    T <- unname(distance(mfdCmp, mu0, muHat)^2 * n)

    if (pvMethod == 'MC') {
      pv <- pvSim(T, ev[ev > 1e-8], pvMC)
    } else {
      pv <- dr::dr.pvalue(ev[ev > 0], T, chi2approx='bx')[['pval.adj']] # Not very accurate
    }
    res <- c(T=T, pv=pv)
  } else if (method == 'proj') {
    if (missing(K)) {
      if (missing(FVE)) {
        stop('K and FVE cannot both be missing')
      } else {
        K <- vapply(FVE, function(fve) {
                      min(which(cumsum(ev) / sum(ev) >= fve))
                       }, 1)
      }
    }
    mfdL2 <- structure(1, class='L2')
    res <- vapply(K, function(k) {
      projs <- sqrt(n) * 
        metric(mfdL2, U=tau(muHat) - tau(mu0), V=ef[, seq_len(k), drop=FALSE])
      T <- unname(sum(projs^2 / ev[seq_len(k)]))
      pv <- stats::pchisq(T, k, lower.tail=FALSE)
      c(T=T, pv=pv)
    }, c(1, 1))
  }

  res
}


pvSim <- function(T, wts, pvMC=1e5) {
  k <- length(wts)
  a <- matrix((stats::rnorm(pvMC * k) ^ 2) * wts, nrow=k)
  pv <- mean(colSums(a) >= T)
  pv
}


# The first partial derivative of the squared geodesic distance. This is defined on the original manifold (without a chart representation)
# V The first argument to D_2 rho^2
# mu The second argument to D_2 rho^2
# Returns: A matrix, with each of the column representing the derivative
D2rho2 <- function(V, mu, project=TRUE) {

  if (!is.matrix(V)) {
    V <- matrix(V)
  }
  ns <- nrow(V)
  xmu <- crossprod(V, mu) / (ns - 1) # Inner product of V and mu
  res <- matrix(-2 * acos(xmu) / sqrt(1 - xmu^2), 
                nrow(V), ncol(V), byrow=TRUE) * 
    V
  if (project) {
    mfdHS <- structure(1, class='HS')
    res <- projectTangent(mfdHS, mu, res)
  }
  res
}


# The second partial derivative of the geodesic distance. This is defined on the original manifold (without a chart representation)
# v The first argument to D_2^2 rho^2
# mu The second argument to D_2^2 rho^2
# Returns: A matrix representing the operator of the bilinear form, operating on the tangent space around the second argument. 
D22rho2 <- function(v, mu, project=TRUE, tol=1e-8) {

  ns <- length(v)
  xmu <- sum(v * mu) / (ns - 1) # Inner product of v and mu
  if (xmu < 1 - tol) {
  res <- -1 / (1 - xmu^2)^(3/2) * 2 * (
    acos(xmu) * xmu * (-1 + xmu^2) * diag(ns) + 
    (acos(xmu) * xmu - sqrt(1 - xmu^2)) * tcrossprod(v) / (ns - 1)
  )
  } else if (xmu >= 1 - tol) {
    res <- 2 * diag(ns)
  }
  if (project) {
    mfdHS <- structure(1, class='HS')
    res <- projectCovS(mfdHS, res, mu)
  }
  res
}


# Get the differential of the inverse chart transition map tau: S -> T_{base}S, taken at tau(mu)
# Push from log_{base}S to log_{mu}S
# This is itself an operator
GetDtauInv <- function(base, mu, mfd, mul=1e-10) {

  stopifnot(inherits(mfd, 'HS'))
  ns <- length(base)
  base <- matrix(base)
  mu <- matrix(mu)
  tau <- function(X) rieLog(mfd, base, X)
  mutau <- tau(mu)

  if (inherits(mfd, 'HS')) {
    if (norm(mfd, U=mutau) <= mul) {
      Dtau <- diag(ns)
    } else {
      normmutau <- norm(mfd, U=mutau)
      Dtau <- sin(normmutau) / normmutau * diag(ns) + 
        normmutau^(-3) * 
          (-normmutau^2 * sin(normmutau) * tcrossprod(base, mutau) / (ns - 1) + 
           (cos(normmutau) * normmutau - sin(normmutau)) * tcrossprod(mutau) / (ns - 1))

    }
  }

  projectTangent(mfd, mu, projMatOnly=TRUE) %*% 
    Dtau %*% 
    projectTangent(mfd, base, projMatOnly=TRUE)
}


# Express the asymptotic covariance function of \hat{mu} for true mean = mu on tau=log_{base}
# First chartless calculation at true mean mu, on the logarithm map based at base:
# Cov(\hmu(s), \hmu(s')) = \Lambda^{-1} CovPsi(X_{i\nu}) \Lambda^{-1}, where
# \Lambda = sample mean of (D_2^2 \rho(X_i, \hmu)^2) as a matrix
# CovPsi(X_{i\nu}) = sample mean of Psi(X_{i\nu})
# This function is implemented on the chart specified by log_base
GetAsympCov <- function(samp, mu, base, mfd) {
  
  nu <- function(X) rieLog(mfd, base, X)
  n <- length(samp$Ly)
  ns <- length(samp$Ly[[1]])
  if (all(sapply(samp$Lt, length) == 1)) {

    X <- matrix(simplify2array(samp$Ly), ncol=n) # Each col is a sample
    PsiX <- GetPsitauX(X, mu, base, mfd)
    # CovPsi <- cov(t(PsiX))

    # Differential of chart transition map
    DtauInv <- GetDtauInv(base, mu, mfd)
    Lambda <- matrix(rowMeans(apply(X, 2, D22rho2, mu=mu, project=FALSE)),
                     ns,
                     ns)
    Lambda <- projectCovS(mfd, Lambda, mu)
    Lambda <- t(DtauInv) %*% Lambda %*% DtauInv
    LamInv <- MASS::ginv(Lambda)
    Z <- LamInv %*% PsiX
    # res <- LamInv %*% CovPsi %*% t(LamInv)
    res <- stats::cov(t(Z))
    res <- (res + t(res)) / 2
    res <- projectCovS(mfd, res, base)

  } else { # longitudinal case for the future
  }

  res
}


# tau = log_{mu0}
# Psitau(X, mu) := D_2 rho_tau^2(X, mu)
GetPsitauX <- function(X, mu, mu0, mfd) {

  t(GetDtauInv(mu0, mu, mfd)) %*% D2rho2(X, mu)


}
