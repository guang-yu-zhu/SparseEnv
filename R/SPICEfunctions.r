## These functions require that the R package
## glasso is installed and loaded


## The following function computes the
## SPICE/(Standardized glasso) precision
## matrix estimate
##
## Arguments:
##   S, the realized sample covariance
##      matrix
##   lam, the positive tuning parameter
## Returns the precision matrix estimate
r.glasso.std=function(S,lam)
{
  p=dim(S)[1]
  std.dev.inv = sqrt(1/diag(S))
  R=std.dev.inv * S * rep(std.dev.inv, each=p);  
  W=glasso::glasso(s=R, rho=lam, penalize.diagonal=FALSE)$wi
  W = std.dev.inv*W*rep(std.dev.inv, each=p)
  return(W)
}


## The following function does cross validation
## to select the tuning parameter for the 
## SPICE/glasso precision matrix estimate
## This function also computes the estimate at 
## the selected tuning parameter.
##
## Arguments:
##   x, the n row by p column matrix where the rows
##      are a realization of n independent copies of
##      a p-variate random vector.
##   lam.vec, a vector of candidate tuning parameter values over which
##            the cross validation procedure searches
##   standard, (TRUE/FALSE), TRUE computes the inverse correlation matrix
##             estimate and then rescales by sample standard deviations
##   ind, optional permutation of 1:n, e.g. ind=sample(n).
##   kfold, the number of folds to use
##
##  Returns a list with the following objects:
##    omega, the inverse covariance estimate computed
##           at the selected tuning parameter value
##    sigma, the inverse of omega
##    best.lam, the selected tuning parameter value
##    cv.err, the vector with jth element equal to the validation error
##            totaled over the kfold folds when the tuning parameter was
##            set to the jth entry in lam.vec (j=1,..., length(lam.vec))
glasso.cv=function(x, lam.vec, standard=TRUE, ind=NULL, kfold=5, quiet = TRUE)
{
  n=dim(x)[1]
  p=dim(x)[2]
  if(is.null(ind)) ind=sample(n);
  cv.loss = array(0, c(length(lam.vec), kfold))
  for (k in 1:kfold)
  {
    leave.out=ind[ (1+floor((k-1)*n/kfold)):floor(k*n/kfold) ]
    x.tr=x[-leave.out,,drop=FALSE]
    meanx=apply(x.tr, 2, mean)
    x.tr=scale(x.tr, center=meanx, scale=FALSE)
    x.va=x[leave.out,,drop=FALSE]
    x.va=scale(x.va, center=meanx, scale=FALSE)
    s.tr=crossprod(x.tr)/(dim(x.tr)[1])
    s.va=crossprod(x.va)/(dim(x.va)[1])
    if(standard)
    {
      std.dev.inv = sqrt(1/diag(s.tr))
      s.tr=std.dev.inv * s.tr * rep(std.dev.inv, each=p);  
    }   
    for(i in 1:length(lam.vec))
    {
      lam=lam.vec[i]
      out=glasso::glasso(s=s.tr, rho=lam, penalize.diagonal=FALSE)
      if(standard)
        omega = std.dev.inv*out$wi*rep(std.dev.inv, each=p) else omega=out$wi;      
      cv.loss[i,k] = sum(omega*s.va) - determinant(omega, logarithm=TRUE)$modulus[1] 
      if(!quiet) cat("Finished lam =", lam.vec[i], "in fold", k, "\n") 
    }
    if(!quiet) cat("Finished fold", k, "\n")       
  }  
  cv.err=apply(cv.loss, 1, sum)
  best.lam = lam.vec[which.min(cv.err)]
  ## compute final estimate at the optimal tuning parameter.
  samp.cov=cov(x)*((n-1)/n)  
  if(standard)
  {
    std.dev.inv = sqrt(1/diag(samp.cov))
    std.dev = sqrt(diag(samp.cov))
    samp.cov=std.dev.inv * samp.cov * rep(std.dev.inv, each=p);  
  }
  out=glasso::glasso(s=samp.cov, rho=best.lam, penalize.diagonal=FALSE)
  if(standard)
  {
    out$wi=std.dev.inv*out$wi*rep(std.dev.inv, each=p)
    out$w=std.dev*out$w*rep(std.dev, each=p)
  }
  return(list(sigma=out$w, omega=out$wi, best.lam=best.lam, cv.err=cv.err)) 
}
