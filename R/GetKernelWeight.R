getKernelWeight <- function(kernel_type,bw,xin,xout,win)
{
    
    invSqrt2pi <-  1/(sqrt(2*pi))
    
    nXGrid <- length(xin);
    nUnknownPoints <- length(xout);
    
    # ========================
    # The checks start here:
    
    if(length(win)!=nXGrid) stop("the length of weight must be same as sample size.")
    
    if(nXGrid == 0) {
        stop("The input X-grid has length zero.")
    }
    
    # Check that bandwidth is greater than zero
    if( bw <= 0){
        stop("The bandwidth supplied for 1-D smoothing is not positive.");
    }
    
    possibleKernels <- c("epan","rect","gauss","gausvar","quar")
    
    # If the kernel_type does not exist, set to epan by default
    if(!(kernel_type %in% possibleKernels)){
        # otherwise use "epan"as the kernel_type 
        warning("Kernel_type argument was not set correctly; Epanechnikov kernel used.");
        Kernel_type = "epan"
    }
    
    
    # Check if the elements are sorted 
    if (is.unsorted(xout)){
        stop("The X-grid used is not sorted.");
    }
    
    # The checks end here.
    # ===================
    kw <- sapply(1:nUnknownPoints,function(i){
        # nested function begins
        index <- (xout[i]-xin<=bw) & (xout[i]-xin>=-bw) # index of elements in the bw window
        if(all(index == FALSE))
            stop('bandwidth is too small for xout[',i,']')
        
        lw <- win[index]
        lx <- xin[index]
        llx = (lx-xout[i]) / bw 
        
        nzw <- switch(kernel_type,
                   epan = (1-llx^2) * 0.75 * lw,
                   rect = lw,
                   gauss = (exp(-0.5*(llx^2))) * invSqrt2pi * lw,
                   gaussvar = lw * (exp(-0.5*llx^2) * invSqrt2pi) *
                       (1.25 - (0.25 * (llx^2))),
                   quar = ((1-llx^2)^2) * (15./16.))
        weight <- rep(0,nXGrid)
        weight[index] <- nzw
        return(weight)
    # nested function end
    })       
    
    return(kw)
}

getKernelWeight1 <- function(kernel_type,bw,xin,xout,win, npoly) {
  n<- length(xin)
  I <- diag(n)  
  kw <- t(sapply(seq_len(n), function(j) {
    Lwls1D(bw, kernel_type, win, xin, I[, j], xout, npoly)
  }))
    
  return(kw)
}
