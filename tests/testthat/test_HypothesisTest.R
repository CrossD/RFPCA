# devtools::load_all()
library(testthat)

test_that('Two sample bootstrap test works for the scalar case', {

            skip('slow')
  nEach <- 25
  n1 <- n2 <- nEach
  mudiff <- 0
  m <- 1
  p <- 40
  tpts <- seq(0, 1, length.out=m)
  spts <- seq(0, 1, length.out=p)
  mfdHS <- structure(1, class=c('HS', 'L2')) # Square-root densities
  mfdL2 <- structure(1, class=c('L2')) # Square-root densities in L2
  mfdDens <- structure(1, class=c('Dens', 'L2')) # Densities, as elements in L2
  sparsity <- p
  K <- 15
  lambda <- 0.05 ^ (seq_len(K) / 3)
  # sum(lambda)

  for (mudiff in c(0, 10)) {
    # mu is a beta pdf 
    a <- 1 + mudiff
    b <- 1
    mu1 <- matrix(dbeta(spts, a, b)) + 1
    mu1 <- sqrt(mu1 / mean(mu1) * (p - 1) / p)
    mu2 <- matrix(dbeta(spts, b, a) + 1)
    mu2 <- sqrt(mu2 / mean(mu2) * (p - 1) / p)
    distance.HS(mfdHS, mu1, mu2)

    set.seed(3)
    samp1 <- MakeMfdProcess(mfdHS, n1, mu1, tpts, K = K, lambda=lambda)
    samp2 <- MakeMfdProcess(mfdHS, n2, mu2, tpts, K = K, lambda=lambda)

    # # Distance from mu
    # hist(distance(mfdHS, matrix(mu1, p, n1), t(drop(samp1$X))))
    # hist(distance(mfdHS, matrix(mu2, p, n2), t(drop(samp2$X))))

    # # Samples are on HS
    # expect_equal(unname(apply(samp1$X, c(1, 3), function(x) sum(x^2) / (p - 1))), matrix(1, n1, m))

    spSamp1 <- SparsifyM(samp1$X, samp1$T, m)
    spSamp2 <- SparsifyM(samp2$X, samp2$T, m)
    yList <- c(spSamp1$Ly, spSamp2$Ly)
    y2List <- sapply(yList, `^`, 2, simplify=FALSE)
    tList <- c(spSamp1$Lt, spSamp2$Lt)
    spSamp <- list(Ly = yList, Lt = tList)
    spSampDens <- list(Ly = y2List, Lt = tList)

    optnsHS <- list(mfd=mfdHS, dataType='Dense', verbose=FALSE)
    optnsL2 <- list(mfd=mfdL2, dataType='Dense', verbose=FALSE)
    optnsDens <- list(mfd=mfdDens, dataType='Dense', verbose=FALSE)

    # resHS <- RFPCA(yList, tList, optnsHS)
    # resDensL2 <- RFPCA(y2List, tList, optnsDens)

    set.seed(1)
    resAll <- BootTest2(spSamp, rep(c(1, 2), each=nEach), mu0=mu1, mfdCmp=mfdHS, RFPCAoptns=optnsHS, parallel='multicore', ncpus=3, method='all', K=c(3, 4))
    set.seed(1)
    res <- BootTest2(spSamp, rep(c(1, 2), each=nEach), mu0=mu1, mfdCmp=mfdHS, RFPCAoptns=optnsHS, parallel='multicore', ncpus=3)
    set.seed(1)
    resProj <- BootTest2(spSamp, rep(c(1, 2), each=nEach), mu0=mu1, mfdCmp=mfdHS, method='proj', K=c(3, 4), RFPCAoptns=optnsHS, parallel='multicore', ncpus=3)
    resPaired <- BootTest2(spSamp, rep(c(1, 2), each=nEach), mu0=mu1, mfdCmp=mfdHS, RFPCAoptns=optnsHS, parallel='multicore', ncpus=3, paired=TRUE)
    resL2 <- BootTest2(spSamp, rep(c(1, 2), each=nEach), mu0=mu1, mfdCmp=mfdHS, RFPCAoptns=optnsL2, parallel='multicore', ncpus=3)
    resDens <- BootTest2(spSampDens, rep(c(1, 2), each=nEach), mu0=mu1, mfdCmp=mfdDens, RFPCAoptns=optnsDens, parallel='multicore', ncpus=3)

    expect_equal(unname(resAll['pv', ]), unname(c(res['pv', ], resProj['pv', ])))
    if (mudiff == 0) {
      expect_true(all(res['pv', ] >= 0.05))
      expect_true(all(resProj['pv', ] >= 0.05))
      expect_true(all(resL2['pv', ] >= 0.05))
      expect_true(all(resDens['pv', ] >= 0.05))
      expect_true(all(resPaired['pv', ] >= 0.05))
    } else if (mudiff >= 1) {
      expect_true(all(res['pv', ] < 0.05))
      expect_true(all(resProj['pv', ] < 0.05))
      expect_true(all(resL2['pv', ] < 0.05))
      expect_true(all(resDens['pv', ] < 0.05))
      expect_true(all(resPaired['pv', ] < 0.05))
    }
  }
})


test_that('One sample bootstrap test works for the scalar case', {

            skip('slow')

  n <- 50
  mudiff <- 0
  m <- 1
  p <- 40
  tpts <- seq(0, 1, length.out=m)
  spts <- seq(0, 1, length.out=p)
  mfdHS <- structure(1, class=c('HS', 'L2')) # Square-root densities
  mfdL2 <- structure(1, class=c('L2')) # Square-root densities in L2
  mfdDens <- structure(1, class=c('Dens', 'L2')) # Densities, as elements in L2
  sparsity <- p
  K <- 15
  lambda <- 0.05 ^ (seq_len(K) / 3)
  # sum(lambda)

  for (mudiff in c(0, 1)) {
    # mu is a beta pdf 
    a <- 1 + mudiff
    b <- 1
    mu0 <- matrix(dbeta(spts, a, b)) + 1
    mu0 <- sqrt(mu0 / mean(mu0) * (p - 1) / p)
    mu <- matrix(dbeta(spts, b, a) + 1)
    mu <- sqrt(mu / mean(mu) * (p - 1) / p)
    distance.HS(mfdHS, mu0, mu)

    set.seed(1)
    samp <- MakeMfdProcess(mfdHS, n, mu, tpts, K = K, lambda=lambda)

    spSamp <- SparsifyM(samp$X, samp$T, m)
    yList <- spSamp$Ly
    y2List <- sapply(yList, `^`, 2, simplify=FALSE)
    tList <- spSamp$Lt
    spSampDens <- list(Ly = y2List, Lt = tList)

    optnsHS <- list(mfd=mfdHS, dataType='Dense', verbose=FALSE)
    optnsL2 <- list(mfd=mfdL2, dataType='Dense', verbose=FALSE)
    optnsDens <- list(mfd=mfdDens, dataType='Dense', verbose=FALSE)

    # resHS <- RFPCA(yList, tList, optnsHS)
    # resDens <- RFPCA(y2List, tList, optnsDens)
    # resDens$muObs

    # Bootstrap test
    set.seed(1)
    resAll <- BootTest1(spSamp, mu0, mfdCmp=mfdHS, RFPCAoptns=optnsHS, parallel='multicore', ncpus=3, method='all', K=c(3, 4))
    set.seed(1)
    res <- BootTest1(spSamp, mu0, mfdCmp=mfdHS, RFPCAoptns=optnsHS, parallel='multicore', ncpus=3)
    # print(system.time({
    set.seed(1)
    resProj <- BootTest1(spSamp, mu0, mfdCmp=mfdHS, method='proj', K=c(3, 4), RFPCAoptns=optnsHS, parallel='multicore', ncpus=3)
    resProjFVE <- BootTest1(spSamp, mu0, mfdCmp=mfdHS, method='proj', FVE=0.99, RFPCAoptns=optnsHS, parallel='multicore', ncpus=3)
    # }))
    resL2 <- BootTest1(spSamp, mu0, mfdCmp=mfdHS, RFPCAoptns=optnsL2, parallel='multicore', ncpus=3)
    resProjL2 <- BootTest1(spSamp, mu0, mfdCmp=mfdHS, method='proj', K=3, RFPCAoptns=optnsL2, parallel='multicore', ncpus=3)

    resDens <- BootTest1(spSampDens, mu0^2, mfdCmp=mfdDens, RFPCAoptns=optnsDens, parallel='multicore', ncpus=3)

    expect_equal(unname(resAll['pv', ]), unname(c(res['pv', ], resProj['pv', ])))
    if (mudiff == 0) {
      expect_true(all(res['pv', ] >= 0.05))
      expect_true(all(resProj['pv', ] >= 0.05))
      expect_true(all(resProjFVE['pv', ] >= 0.05))
      expect_true(all(resL2['pv', ] >= 0.05))
      expect_true(all(resProjL2['pv', ] >= 0.05))
    } else if (mudiff > 1 - 1e-5) {
      expect_true(all(res['pv', ] < 0.05))
      expect_true(all(resProj['pv', ] < 0.05))
      expect_true(all(resProjFVE['pv', ] < 0.05))
      expect_true(all(resL2['pv', ] < 0.05))
      expect_true(all(resProjL2['pv', ] < 0.05))
      expect_true(all(resDens['pv', ] < 0.05))
    }
  }
})


test_that('Helper functions for the asymptotic test works', {

  n <- 50 # Test small n
  mudiff <- 1
  m <- 1
  p <- 400
  tpts <- seq(0, 1, length.out=m)
  spts <- seq(0, 1, length.out=p)
  mfdHS <- structure(1, class=c('HS', 'L2')) # Square-root densities
  sparsity <- p
  K <- 15
  lambda <- 0.05 ^ (seq_len(K) / 3)

  for (mudiff in c(0, 1)) {
    # mu is a beta pdf 
    a <- 1 + mudiff
    b <- 1
    mu0 <- matrix(dbeta(spts, a, b)) + 1
    mu0 <- sqrt(mu0 / mean(mu0) * (p - 1) / p)
    mu <- matrix(dbeta(spts, b, a) + 1)
    mu <- sqrt(mu / mean(mu) * (p - 1) / p)
    manifold:::distance.HS(mfdHS, mu0, mu)

    set.seed(1)
    samp <- MakeMfdProcess(mfdHS, n, mu, tpts, K = K, lambda=lambda)

    spSamp <- SparsifyM(samp$X, samp$T, m)
    yList <- spSamp$Ly
    tList <- spSamp$Lt
    optnsHS <- list(mfd=mfdHS, dataType='Dense', verbose=FALSE)

    # GetDtauInv: Change chart from log_mu0 to log_mu
    Dtau <- GetDtauInv(mu0, mu, mfd=mfdHS)
    delta <- 1e-8
    x <- mu
    y <- samp$X[1, , ]
    v <- rieLog(mfdHS, x, y)
    xdelta <- rieExp(mfdHS, x, v * delta)
    tauFrom <- function(X) rieLog(mfdHS, mu0, X)
    tauTo <- function(X) rieLog(mfdHS, mu, X)
    fromCoord <- tauFrom(cbind(x, xdelta))
    toCoord <- tauTo(cbind(x, xdelta))
    hFrom <- matrix((fromCoord[, 2] - fromCoord[, 1]) / delta)
    hTo <- matrix((toCoord[, 2] - toCoord[, 1]) / delta)
    trans <- Dtau %*% hFrom
    # Dtau successfully change the coordinates 
    expect_equal(trans, hTo, tolerance=1e-6)
    expect_equal(c(crossprod(trans, mu)), 0, tolerance=1e-5)
    expect_equal(max(abs(crossprod(Dtau, mu))), 0)

    if (mudiff == 0) {
      ev <- eigen(Dtau)$values
      ev[p] <- ev[p] + 1
      expect_equal(ev, rep(1, p)) # All but the last eigenvalue of Dtau is 1
    } else if (mudiff >= 1) {
      expect_equal(max(abs(
                     t(mu) %*% Dtau %*% rieLog(mfdHS, mu0, t(drop(samp$X)))
                   )), 0)
    }

    # D2rho2
    dProj1 <- D2rho2(samp$X[1, , ], mu, project=TRUE)
    dProjAll <- D2rho2(t(drop(samp$X)), mu, project=TRUE)
    expect_equal(c(dProjAll[, 1]), c(dProj1))
    expect_equal(max(abs(crossprod(dProjAll, mu))), 0)


    # D22rho2
    x <- samp$X[1, , ]
    v <- rieLog(mfdHS, mu, x)
    v <- Normalize(v) * sqrt(p - 1)
    d2ProjOrth <- D22rho2(v, mu, project=TRUE)
    ev <- eigen(d2ProjOrth)$values
    expect_equal(ev, c(2, rep(0, p - 1)))

    d2Projmu <- D22rho2(mu, mu, project=TRUE)
    ev1 <- eigen(d2Projmu)$values
    expect_equal(ev1, c(rep(2, p - 1), 0))

    d2Proj <- D22rho2(x, mu, project=TRUE)
    d2NoProj <- D22rho2(samp$X[1, , ], mu, project=FALSE)
    V <- rieLog(mfdHS, mu, t(drop(samp$X)))
    valProj <- t(V) %*% d2Proj %*% V / p^2
    valNoProj <- t(V) %*% d2NoProj %*% V / p^2
    # As bilinear operators operating on tangent vectors, projection does nothing
    expect_equal(valProj, valNoProj)
    # As an operator, projection is needed
    vProj <- d2Proj %*% V / p
    vNoProj <- d2NoProj %*% V / p
    # The intrinsic bilinear operator is on the tangent space at mu
    expect_equal(max(abs(crossprod(vProj, mu))), 0)
    expect_true(max(abs(crossprod(vNoProj, mu))) > 1e-10)

    acov <- GetAsympCov(spSamp, mu, mu0, mfdHS)
    expect_equal(max(abs(acov %*% mu0)), 0)


  }
})


test_that('Two sample asymptotic test works for the scalar case', {

            skip('slow')

  MC <- 50
  nEach <- 50
  n1 <- n2 <- nEach
  mudiff <- 0
  m <- 1
  p <- 40
  tpts <- seq(0, 1, length.out=m)
  spts <- seq(0, 1, length.out=p)
  mfdHS <- structure(1, class=c('HS', 'L2')) # Square-root densities
  sparsity <- p
  K <- 15
  lambda <- 0.05 ^ (seq_len(K) / 3)
  # sum(lambda)

  for (mudiff in c(0, 10)) {
    # mu is a beta pdf 
    a <- 1 + mudiff
    b <- 1
    mu1 <- matrix(dbeta(spts, a, b)) + 1
    mu1 <- sqrt(mu1 / mean(mu1) * (p - 1) / p)
    mu2 <- matrix(dbeta(spts, b, a) + 1)
    mu2 <- sqrt(mu2 / mean(mu2) * (p - 1) / p)
    distance.HS(mfdHS, mu1, mu2)

    optnsHS <- list(mfd=mfdHS, dataType='Dense', verbose=FALSE)

    # resPaired <- AsympTest2(spSamp, rep(c(1, 2), each=nEach), mu0=mu1, mfdCmp=mfdHS, paired=TRUE, RFPCAoptns=optnsHS)

    res <- sapply(seq_len(MC), function(i) {
      set.seed(i)
      samp1 <- MakeMfdProcess(mfdHS, n1, mu1, tpts, K = K, lambda=lambda)
      samp2 <- MakeMfdProcess(mfdHS, n2, mu2, tpts, K = K, lambda=lambda)
      spSamp1 <- SparsifyM(samp1$X, samp1$T, m)
      spSamp2 <- SparsifyM(samp2$X, samp2$T, m)
      yList <- c(spSamp1$Ly, spSamp2$Ly)
      # y2List <- sapply(yList, `^`, 2, simplify=FALSE)
      tList <- c(spSamp1$Lt, spSamp2$Lt)
      spSamp <- list(Ly = yList, Lt = tList)
      # spSampDens <- list(Ly = y2List, Lt = tList)

      resAsymp <- AsympTest2(spSamp, rep(c(1, 2), each=nEach), mu0=mu1, mfdCmp=mfdHS, RFPCAoptns=optnsHS)
      resAsympProj <- AsympTest2(spSamp, rep(c(1, 2), each=nEach), mu0=mu1, mfdCmp=mfdHS, RFPCAoptns=optnsHS, method='proj', K=3)
      cbind(resAsymp, resAsympProj)
    }, simplify='array')
    pv <- res['pv', , ]
    # hist(pv)

    if (mudiff == 0) {
      expect_true(all(rowMeans(pv <= 0.05) <= 0.1))
    } else if (mudiff >= 1) {
      expect_true(all(rowMeans(pv <= 0.05) > 0.6))
    }
  }
})


test_that('One sample asymptotic test works for the scalar case', {

            skip('slow')

  MC <- 50
  n <- 50
  mudiff <- 0
  m <- 1
  p <- 40
  tpts <- seq(0, 1, length.out=m)
  spts <- seq(0, 1, length.out=p)
  mfdHS <- structure(1, class=c('HS', 'L2')) # Square-root densities
  mfdL2 <- structure(1, class=c('L2')) # Square-root densities in L2
  mfdDens <- structure(1, class=c('Dens', 'L2')) # Densities, as elements in L2
  sparsity <- p
  K <- 15
  lambda <- 0.05 ^ (seq_len(K) / 3)
  # sum(lambda)

  for (mudiff in c(0, 1)) {
    # mu is a beta pdf 
    a <- 1 + mudiff
    b <- 1
    mu0 <- matrix(dbeta(spts, a, b)) + 1
    mu0 <- sqrt(mu0 / mean(mu0) * (p - 1) / p)
    mu <- matrix(dbeta(spts, b, a) + 1)
    mu <- sqrt(mu / mean(mu) * (p - 1) / p)
    distance.HS(mfdHS, mu0, mu)

    i <- 1
    res <- sapply(seq_len(MC), function(i) {
      set.seed(i)
      samp <- MakeMfdProcess(mfdHS, n, mu, tpts, K = K, lambda=lambda)

      spSamp <- SparsifyM(samp$X, samp$T, m)
      yList <- spSamp$Ly
      tList <- spSamp$Lt

      optnsHS <- list(mfd=mfdHS, dataType='Dense', verbose=FALSE)

      # Asymptotic test
      # L2 distance
      resAsymp <- AsympTest1(spSamp, mu0, mfdCmp=mfdHS, RFPCAoptns=optnsHS, pvMethod='MC')
      # Projection
      resAsymp1 <- AsympTest1(spSamp, mu0, mfdCmp=mfdHS, method='proj', K=3, RFPCAoptns=optnsHS, pvMethod='MC')
      resAsymp2 <- AsympTest1(spSamp, mu0, mfdCmp=mfdHS, method='proj', FVE=0.9, RFPCAoptns=optnsHS, pvMethod='MC')
      cbind(resAsymp, resAsymp1, resAsymp2)
    }, simplify='array')
    pv <- res['pv', , ]

    if (mudiff == 0) {
      expect_true(all(rowMeans(pv <= 0.05) <= 0.1))
    } else if (mudiff > 1 - 1e-5) {
      expect_true(all(rowMeans(pv <= 0.05) >= 0.7))
    }
  }
})

