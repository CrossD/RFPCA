#' @export
metric <- function(mfd,U,V) UseMethod("metric",mfd)

#' @export
norm <- function(mfd,U) UseMethod("norm",mfd)

#' @export
distance <- function(mfd,X,Y) UseMethod("distance",mfd)

#' @export
geodesicCurve <- function(mfd,p,h,t) UseMethod("geodesicCurve",mfd)

#' @export
rieExp <- function(mfd,p,V) UseMethod("exp",mfd)

#' @export
rieLog <- function(mfd,p,X) UseMethod("log",mfd)

#' @export
# Project a matrix whose columns are data points onto the manifold
project <- function(mfd,p) UseMethod("project",mfd)

geopolate <- function(mfd,...) UseMethod("geopolate",mfd)

#' @export
projectTangent <- function(mfd, ...) UseMethod('projectTangent', mfd)

rmfd <- function(mfd,n,...) UseMethod("rmfd",mfd)

# dotprod <- function(u,v) {
    # return(colSums(u*v))
# }

# dimAmbient Dimension of the ambient space.
# dimIntrinsic Intrinsic dimension of the manifold.
# dimTangent Number of parameters for the tangent space, or the length of a tangent vector
calcGeomPar <- function(mfd, dimAmbient, dimIntrinsic, dimTangent) {
  if (sum(c(!missing(dimIntrinsic), !missing(dimAmbient), !missing(dimTangent))) > 1) {
    stop('Can only specify one of dimIntrinsic, dimAmbient, and dimTangent')
  }
  UseMethod('calcGeomPar', mfd)
}
calcIntDim <- function(mfd, dimAmbient, dimTangent) {
  if (!missing(dimAmbient) && !missing(dimTangent)) {
    stop('Cannot specify both dimAmbient and dimTangent')
  }
  UseMethod('calcIntDim', mfd)
}
calcTanDim <- function(mfd, dimAmbient, dimIntrinsic) {
  if (!missing(dimAmbient) && !missing(dimIntrinsic)) {
    stop('Cannot specify both dimAmbient and dimIntrinsic')
  }
  UseMethod('calcTanDim', mfd)
}


frechetMean <- function(mfd, X, weight=NULL, tol=1e-9, maxit=1000, lam=0.99) {

  X <- as.matrix(X)
  dimAmbient <- nrow(X)
  dimTangent <- calcTanDim(mfd, dimAmbient)
  
  if (is.null(weight)) {
    n <- dim(X)[2]
    weight <- rep(1 / n, n)
    mu0 <- apply(X, 1, mean)
  } else {
    weight <- weight / sum(weight)
    ind <- weight > 1e-14
    n <- sum(ind)
    X <- X[, ind, drop=FALSE]
    weight <- weight[ind]
    mu0 <- apply(X*t(replicate(dimAmbient, weight)), 1, sum)
  }

  it <- 0
  dif <- Inf
  SS <- 0

  # browser()
  while(dif > tol && it < maxit) {
    mu0 <- project(mfd, mu0)
    it <- it + 1
    V <- rieLog(mfd, mu0, X)
    vNew <- as.numeric(rowSums(V * matrix(weight, dimTangent, n, byrow=TRUE)))
    muNew <- as.numeric(rieExp(mfd, mu0, vNew * lam^it))
    SSNew <- sum(V^2)
    dif <- abs(SS - SSNew) / (SS + tol)

    mu0 <- muNew
    SS <- SSNew
  # browser()
  }

  if (it == maxit) {
    stop('Maximum iteration reached!')
  }

  matrix(mu0)
}



# #' Frechet mean curve of a set of sparsely observed curves on the manifold
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param bw The bandwidth
# #' @param kernel_type The type of kernel for smoothing
# #' @param yin A list of \emph{n} matrices containing the observed values for each individual. Missing values specified by \code{NA}s are supported for dense case (\code{dataType='dense'}).
# #' @param xin A list of \emph{n} vectors containing the observation time points for each individual corresponding to y. Each vector should be sorted in ascending order.
# #' @param xout The time points where the mean curve shall be evaluated
# #' @param npoly The order of the local polynomial smoother
# #' @param win A vector \emph{n} numbers. The weight of each individual.
# #' 
# #' @return A matrix with \code{length(xout)} columns. The \emph{i}th column is the estimated mean curve at \code{xout[i]} 

frechetMeanCurve <- function(mfd, bw, kernel_type, xin, 
                             yin, xout, npoly = 1L, 
                             win=rep(1L, length(xin))){
  
  if(is.unsorted(xout)){
      stop('`xout` needs to be sorted in increasing order')
  }
  
  if(all(is.na(win)) || all(is.na(xin)) || all(is.na(yin))){
      stop(' win, xin or yin contain only NAs!')
  }
  
  # Deal with NA/NaN measurement values
  NAinY = is.na(xin) | is.na(yin) | is.na(win)
  if(any(NAinY)){
      win = win[!NAinY]
      xin = xin[!NAinY]
      yin = yin[!NAinY]
  } 
  
  if (!inherits(mfd, 'Euclidean')) {
    nout <- length(xout)
    # kw <- getKernelWeight(kernel_type,bw,xin,xout,win)
    # browser()
    kw <- getKernelWeight1(kernel_type,bw,xin,xout,win, npoly)
    meanCurve <- sapply(1:nout, function(i) frechetMean(mfd,yin,weight=kw[,i]))
    
    # if(npoly == 0) # constant kernel smoothing
    # {
        # could be faster if local compact kernel by drop points with zero weight
    # }
    # else # local linear smoothing
    # {
        # stop('npoly >= 1 not supported yet')
        # #meanCurve <- matrix(0,3,nout)
    # }
  } else if (inherits(mfd, 'Euclidean')) {
    meanCurve <- t(apply(yin, 1, function(x) {
                           Lwls1D(bw, kernel_type, win, xin, x, xout, npoly)
                   }))
  }
  return(meanCurve)
}


# This function computes the optimal bandwidth choice for the mean
# function use GCV method by pooling the longitudinal data together. 
# verbose is unused for now
# this is compatible with PACE because the GCV is calculated in a same way

GCVFrechetMeanCurve <- function(mfd, yy, tt, kernel, npoly, dataType, verbose=TRUE) {
  
  # If 'yy' and 't' are vectors "cheat" and break them in 
  # a list of 10 elements                                                                                                                                  
  if ( is.vector(yy) && is.vector(tt) && !is.list(yy) && !is.list(tt) ){
    if (length(tt) < 21) {
        stop("You are trying to use a local linear weight smoother in a vector with less than 21 values.\n")
    }
    myPartition = c(1:10, sample(10, length(tt)-10, replace=TRUE));
    yy = split(yy, myPartition)
    tt = split(tt, myPartition)
    dataType = 'Sparse';
  }
  
  # pool all observations across subjects and re-order them so that
  # t is ascending
  t = unlist(tt);
  y = do.call(cbind, lapply(yy, unlist)) # unlist(yy)[order(t)];
  torder <- order(t)
  t = t[torder]
  y = y[, torder]
  
  # r = diff(range(t))
  N = length(t);
  r = t[N] - t[1];  
  
  # Specify the starting bandwidth candidates
  if (dataType == "Sparse") {
      dstar = fdapace:::Minb(t, npoly+2)
      if ( dstar > r*0.25){
          dstar = dstar * 0.75
          warning( c( "The min bandwidth choice is too big, reduce to ", dstar, "!\n"))
      }
      h0 = 2.5 * dstar
  } else if (dataType == "DenseWithMV") {
      h0 = 2.0 * fdapace:::Minb(t, npoly+1)
  } else {
      h0 = 1.5 * fdapace:::Minb(t, npoly+1)
  }  
  if ( is.nan(h0) ){
      if ( kernel == "gauss" ){
          h0 = 0.2 * r
      } else {
          stop("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")
      }
  }  
  h0 <- min(h0, r)
  q = (r / (4 * h0))^(1/9);
  bwCandidates = sort(q^(0:9)*h0)
  
  # idx = apply(X= sapply(X=t, FUN='==',  ...=sort(unique(t)) ),MARGIN=2, FUN=which)  
  idx =  pracma::uniq(t)$n
  # This is to make sure we get the same as MATLAB PACE 
  # I would write them in a function (equivalent of mykernel.m) if it is worth it 
  # Similarly there is no reason to repeat the FOR-loop twice; this too can go into a seperate function
  k0_candidates <- list('quar' = 0.9375,  'epan' = 0.7500, 'rect' = 0.5000, 
                        'gausvar' = 0.498677, 'gausvar1' = 0.598413,  'gausvar2' = 0.298415, 'other' = 0.398942)
  if( any(names(k0_candidates) == kernel)){ 
      k0 = as.numeric(k0_candidates[kernel])
  } else { 
      k0 =  as.numeric(k0_candidates$other)
  }
  
  gcvScores <- c()

  # Get the corresponding GCV scores 
  for(i in 1:length(bwCandidates)) {
      # newmu = Lwls1D(bwCandidates[i], kern=kernel, npoly=npoly, nder=nder, xin = t,yin= y,xout= sort(unique(t)))[idx]
    # browser()
    # if (!inherits(mfd, 'Euclidean')) {
      # newmu = frechetMeanCurve(mfd, bwCandidates[i], 
                               # kernel_type=kernel, 
                               # npoly=npoly, 
                               # xin = t,
                               # yin= y, 
                               # win = rep(1,dim(y)[2]),
                               # xout= sort(unique(t)))[,idx]
    # } else if (inherits(mfd, 'Euclidean')) {
      # newmu <- t(apply(y, 1, function(x) {
          # Lwls1D(bwCandidates[i], kernel_type=kernel, xin=t, yin=x, npoly=npoly, xout=sort(unique(t)))
      # }))[, idx]
    # }
    newmu <- frechetMeanCurve(mfd, bwCandidates[i], 
                              kernel_type=kernel, 
                              npoly=npoly, 
                              xin = t,
                              yin= y, 
                              win = rep(1,dim(y)[2]),
                              xout= sort(unique(t)))[,idx]
    cvsum = sum((newmu -y)^2)
    gcvScores[i] =cvsum/(1-(r*k0)/(N*bwCandidates[i]))^2
  }
  
  # If no bandwith gives a finite gcvScore increase the candidate bandwith and retry on a finer grid
  if (all((is.infinite(gcvScores)))) {
      bwCandidates = seq( max(bwCandidates), r, length.out = 2*length(bwCandidates))
      for(i in 1:length(bwCandidates)){
          # newmu = Lwls1D(bwCandidates[i], kern=kernel, npoly=npoly, nder=nder, xin = t,yin= y,xout= sort(unique(t)))[idx]
        # if (!inherits(mfd, 'Euclidean')) {
        newmu <- frechetMeanCurve(mfd, bwCandidates[i], 
                                  kernel_type=kernel, 
                                  npoly=npoly, 
                                  xin = t,
                                  yin= y, 
                                  win = rep(1,dim(y)[2]),
                                  xout= sort(unique(t)))[,idx]
        # } else if (inherits(mfd, 'Euclidean')) {
          # newmu <- t(apply(y, 1, function(x) {
                             # Lwls1D(bwCandidates[i], kernel_type=kernel, xin=t, yin=x, 
                                    # npoly=npoly, xout=sort(unique(t)))
                           # }))[, idx]
        # }
          cvsum = sum((newmu -y)^2 )
          gcvScores[i] =cvsum/(1-(r*k0)/(N*bwCandidates[i]))^2
      }
  }
  
  # If the problem persist we clearly have too sparse data
  if(all((is.infinite(gcvScores)))){
    stop("The data is too sparse, no suitable bandwidth can be found! Try Gaussian kernel instead!\n")      
  }
  
  bInd = which(gcvScores == min(gcvScores));
  bScr = gcvScores[bInd][1]
  bOpt = max(bwCandidates[bInd]); 
  
  if( bOpt == r){
    warning("data is too sparse, optimal bandwidth includes all the data!You may want to change to Gaussian kernel!\n")
  }
  bOptList <- list( 'bOpt' = bOpt, 'bScore' = bScr) 
  return( bOptList)
}

