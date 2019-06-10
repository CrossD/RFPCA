# #' Riemannian metric of tangent vectors
# #' 
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param p The base point which could be NULL if the metric does not rely on the base points
# #' @param U A D*n matrix, each column represents a tangent vector at \emph{p}
# #' @param V A D*n matrix, each column represents a tangent vector at \emph{p}
# #' 
# #' @details The tangent vectors can be represented in a coordinate frame, or in ambient space
metric.Euclidean <- function(mfd,U,V,p=NULL)
{
  U <- as.matrix(U)
  V <- as.matrix(V)
  return(colSums(U * V))
}

# #' The norm induced by the Riemannian metric tensor on tangent spaces
# #' 
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param p The base point which could be NULL if the norm does not rely on it
# #' @param U A D*n matrix, where each column represents a tangent vector at \emph{p}
norm.Euclidean <- function(mfd,U,p=NULL)
{
    return(sqrt(colSums(U*U)))
}

# #' Geodesic distance of points on the manifold
# #' 
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param X A D*n matrix. Each column represents a point on manifold
# #' @param Y A D*n matrix. Each column represents a point on manifold
# #' @return A 1*n vector. The \emph{i}th element is the geodesic distance of \code{X[,i]} and \code{Y[,i]}
distance.Euclidean <- function(mfd,X,Y)
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  sqrt(colSums((X - Y)^2))
}

# #' Riemannian exponential map at a point p
exp.Euclidean <- function(mfd, p, V) {

  V <- as.matrix(V)

  res <- V + matrix(p, nrow=nrow(V), ncol=ncol(V))
  res
}


# #' Riemannian log map at a point
# #' @param mfd A class instance that represents the Riemannian manifold
# #' @param p A matrix. Each column represents a base point of the log map. If only one base point is supplied, then it is replicated to match the number of points in \emph{X}
# #' @param X A matrix. Each column represents a point on the manifold
# #' @return A matrix with the \emph{i}th column being the log map of the \emph{i}th point

log.Euclidean <- function(mfd, p, X)
{

  X <- as.matrix(X)
  Z <- X - matrix(p, nrow(X), ncol(X))

  return(Z)
}

#' @export
project.Euclidean <- function(mfd, A) A


calcGeomPar.Euclidean <- function(mfd, dimAmbient, dimIntrinsic, dimTangent) {

  if (!missing(dimAmbient)) {
    dimAmbient
  } else if (!missing(dimIntrinsic)) {
    dimIntrinsic
  } else if (!missing(dimTangent)) {
    dimTangent
  }
}


calcIntDim.Euclidean <- function(mfd, dimAmbient, dimTangent) {
  
  if (!missing(dimAmbient)) {
    dimAmbient
  } else if (!missing(dimTangent)) {
    dimTangent
  }
}


calcTanDim.Euclidean <- function(mfd, dimAmbient, dimIntrinsic) {

  if (!missing(dimAmbient)) {
    dimAmbient
  } else if (!missing(dimIntrinsic)) {
    dimIntrinsic
  }
}


# Project ambient space data onto the tangent space at p
projectTangent.Euclidean <- function(mfd, p, X, projMatOnly=FALSE) {

  if (projMatOnly) {
    return(diag(length(p)))
  } else {
    return(X)
  }
}
