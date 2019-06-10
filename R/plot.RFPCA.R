#' Plot an RFPCA object
#'
#' @param x An RFPCA object
#' @param type The type of plot to produce. 'mu', mean function; 'phi', eigenfunctions; 'cov', covariance function'; 'fitted', fitted trajectories; 'raw', raw observations
#' @param K The number of components to use or plot
#' @param subset A subset of subject IDs to plot
#' @param inputData The input data (a list of Ly and Lt)
#' @param dimNames An optional vector specifying the dimension names
#' @export

plot.RFPCA <- function(x, type=c('mu', 'phi', 'cov', 'fitted', 'raw'), K=3, subset, inputData, dimNames) {

  type <- match.arg(type)
  if (!missing(dimNames)) {
    stopifnot(length(dimNames) == nrow(x[['muWork']]))
  } else if (missing(dimNames)) {
    dimNames <- dimnames(x[['muWork']])[[1]]
    if (is.null(dimNames)) {
      dimNames <- paste0('V', seq_len(nrow(x[['muWork']])))
    }
  }
  if (dim(x[['phi']])[3] < K & !type %in% c('mu', 'cov', 'raw')) {
    stop(sprintf('Only %d components are available; cannot print K = %d components', dim(x[['phi']])[3], K))
  }

  tt <- x[['workGrid']]


  if (type == 'mu') {
    dimnames(x[['muWork']])[[1]] <- dimNames
    dimnames(x[['muWork']])[[2]] <- tt
    muDat <- cbind(melt(x[['muWork']], varnames=c('variable', 't')))
    p <- ggplot(muDat, aes(x=t, y=value, group=variable, color=variable)) + geom_line(size=2) + ylab(expression(hat(mu(t))))
  } else if (type == 'phi') {
    phiDat <- structure(x[['phi']], dimnames=
                             list(t = tt, 
                                  variable = dimNames, 
                                  K = seq_len(x[['K']]))) 
    phiDat <- phiDat[, , seq_len(min(K, dim(phiDat)[3])), drop=FALSE]
    phiDat <- melt(phiDat)
    lam <- x[['lam']]
    phiDat$Kname <- sprintf("phi[%s] ~ ~ %.4g * \'%%\'", phiDat$K, lam[phiDat$K] / sum(lam) * 100)
    p <- ggplot(phiDat, aes(x=t, y=value, color=variable)) + geom_line(size=2) + facet_grid( ~ Kname, scales='free_y', labeller=label_parsed) + ylab(expression(hat(phi(t))))
  } else if (type == 'cov') {
    covDat <- structure(x[['cov']], dimnames=
                             list(t1 = tt, 
                                  t2 = tt, 
                                  variable1 = dimNames, 
                                  variable2 = dimNames)) 
    covDat <- melt(covDat)
    p <- ggplot(covDat, aes(x=t1, y=t2, fill=value)) + geom_raster() + facet_grid(variable1 ~ variable2)
  } else if (type == 'fitted') {
    fitPool <- fitted(x, K, grid='work')
    n <- dim(fitPool)[1]
    dimnames(fitPool) <- list(id = seq_len(n), 
                              variable = dimNames, 
                              t = tt)
    fitPool <- melt(fitPool)
    if (!missing(subset)) {
      fitPool <- fitPool[fitPool$id %in% subset, , drop=FALSE]
    } else {
      subset <- seq_len(n)
    }
    p <- ggplot(fitPool, aes(x=t, y=value, color=variable)) + geom_line() + facet_wrap(~id)

    if (!missing(inputData)) {
      tRaw <- inputData$Lt[subset]
      datRaw <- data.frame(do.call(rbind, lapply(inputData$Ly[subset], t)))
      names(datRaw) <- dimNames
      datRaw <- cbind(datRaw, t=unlist(tRaw), id = rep(subset, vapply(tRaw, length, 1L)))
      datRaw <- melt(datRaw, c('t', 'id'))
      p <- p + geom_point(data=datRaw, aes(x=t, y=value, color=variable)) + facet_wrap(~id)
    }

  }

  print(p)
  invisible(p)
}
