#' RFPCA: Riemannian Functional Principal Component Analysis
#'
#' RFPCA for Functional Data Analysis of Riemannian manifold data
#' 
#' References: 
#' Dai X, Lin Z, Müller HG. Modeling sparse longitudinal data on Riemannian manifolds. Biometrics. 2021;77(4):1328–41.
#' Lin Z, Yao F. Intrinsic Riemannian functional data analysis. The Annals of Statistics. 2019;47(6):3533–77.
#' Dai X, Müller HG. Principal component analysis for functional data on Riemannian manifolds and spheres. Annals of Statistics. 2018;46(6B):3334–61.
#' 
#'  
#' 
#' Maintainer:  Xiongtao Dai \email{xdai@@iastate.edu}
#' 
#' @author
#' Xiongtao Dai \email{xdai@@iastate.edu}
#' Zhenhua Lin
#'
#' 
#'
#' @docType package
#' @name RFPCA
#' @import manifold
#' @importFrom stats setNames pchisq prcomp cov rnorm rexp rt fitted kernel
#' @importFrom fdapace FPCA ConvertSupport Lwls1D Lwls2D
#' @importFrom utils methods
#' @importFrom abind abind
NULL
