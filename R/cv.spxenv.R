cv.spxenv <- function (X, Y, dims=NULL,lambda=NULL,fold = 5, maxiter=1e2,ftol=1e-3,seed=2015,paralell=FALSE,core=2) 
{
  time1=proc.time()
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(Y)
  if(is.null(dims)) dims=1:p
  if(missing(lambda)) lambda <- exp(seq(log(1),log(0.001),len=5))
  #foldi<-list()
  #for(i in 1:fold) {
  #  foldi[[i]]=(floor((i - 1) * n / fold) + 1) : ceiling(i * n / fold)
  #}
  set.seed(seed)
  foldi <- split(sample(1:n), rep(1:fold, length = n))
  mspemati = matrix(0,fold,length(dims),dimnames=list(fold=1:fold,num0=dims))
  for (i in 1:fold) {
    omit <- foldi[[i]]
    tempX <- as.matrix(X[-omit,])
    tempY <- as.matrix(Y[-omit,])
    testX <- as.matrix(X[omit,])
    testY <- as.matrix(Y[omit,])
    for(j in 1 : length(dims)){ 
      cat('----------fold=',i,',dims=',dims[j],'\n')
      m0 <- spxenv(tempX,tempY,dims[j],lambda=lambda,maxiter=maxiter,ftol=ftol)
      print(m0$where1)
      pred <- matrix(1, nrow(testX),1) %*% t(m0$mu) + testX %*% m0$beta
      resi <- testY - pred
      mspemati[i,j] <- sqrt(mean(apply(resi^2, 1, sum)))        
    }        
  }
  mspemat<- apply(mspemati, 2, mean)
  u.opt = dims[which.min(mspemat)]  
  cv <- list(mspemat = mspemat,u.opt=u.opt,fit.time=(proc.time()-time1)[3])
  cv
}

