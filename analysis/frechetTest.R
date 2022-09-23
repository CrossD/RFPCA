# Code shared by Dr Paromita Dubey
# Reference: Dubey P, Müller H-G. Fréchet analysis of variance for random objects. Biometrika. 2019 Dec 1;106(4):803–21. 


#Y1: list of covariance matrix inputs of first sample
#Y2: list of covariance matrix inputs of second sample
#fmean(Y): function to compute frechet mean of the elements in the list Y
#Dist(Y,A): function to compute squared distances of the elements in the list Y from the matrix A. returns a 
#vector of the same length as Y


#computes the frechet variance of the elements in A around the matrix in B
frechet_var<-function(A,B) {
  return(mean(Dist(A,B)))
}


#computes the asymptotic variance of the frechet variance of the elements in A around the matrix B
sigma <-function(A,B){
  D=Dist(A,B)
  return(mean(D^2)-(mean(D))^2)
}


frechet.statistic<- function(sizes, data, indices) { 
  L<- data[indices] # allows boot to select sample 
  n1 <-sizes[1]
  n2<-sizes[2]
  n=n1+n2
  lambda=n1/n
  A<-L[1:n1]
  B<-L[(n1+1):n]
  A_mean=fmean(A)
  B_mean=fmean(B)
  L_mean=fmean(L)
  V1=frechet_var(A,A_mean)
  V2=frechet_var(B,B_mean)
  V1n=frechet_var(A,B_mean)
  V2n=frechet_var(B,A_mean)
  V=frechet_var(L,L_mean)
  sigma0=sigma(L,L_mean)
  add.factor<-(sqrt(n)*(V1n-V1))+(sqrt(n)*(V2n-V2))
  Test.asy<- ((lambda*(1-lambda)*((n*(V1-V2)^2)+add.factor^2)/sigma0))
  return(Test.asy)
}


test<-function(Y1,Y2, R=1000, nCores=1L) {
  L=append(Y1,Y2)
  n1=length(Y1)
  n2=length(Y2)
  f.results <- boot::boot(data=L, statistic=frechet.statistic, sim = "ordinary",
                    R=R, sizes=c(n1,n2), parallel='multicore', ncpus = nCores)
  pval.asy=1-pchisq(f.results$t0,1)
  pval.boot= length(which(f.results$t>f.results$t0))/R
  return(list(pval.asy=pval.asy,pval.boot=pval.boot))
}


fmean <- function(Y) {
  Y <- matrix(unlist(Y), ncol=length(Y))
  res <- as.numeric(cov2cor(matrix(frechetMean(mfd, Y), d, d)))
  matrix(res, d, d)
}


Dist <- function(Y, A) {
  Y <- matrix(unlist(Y), ncol=length(Y))
  distance(mfd, matrix(A, d^2, ncol(Y)), Y)
}


