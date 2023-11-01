spxenv_spice<-function(X,Y,lambda_spice=0.1){
  if(is.null(lambda_spice))  lambda_spice=0.1
  spice = list()
  n = nrow(X)
  p = ncol(X)
  r = ncol(Y)
  if(r<n){ 
    Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
    Xc <- as.matrix(scale(X, center = T, scale = FALSE))  
    sigX <- stats::cov(Xc)*(n-1)/n
    sigY <- stats::cov(Yc)*(n-1)/n 
    sigXY <- stats::cov(Xc, Yc)*(n-1)/n
    invsigY <- chol2inv(sechol(sigY)) # invsigX=tmp2%*%t(tmp2)
    sigXcY <- sigX-sigXY%*% invsigY %*% t(sigXY)
    
    invsigXcY_spice= r.glasso.std(sigXcY,lam = lambda_spice)
    spice$invsigXcY <- invsigXcY_spice
    spice$sigXcY <- chol2inv(chol(invsigXcY_spice))
    
    invsigX_spice= r.glasso.std(sigX,lam = lambda_spice)
    spice$invsigX = invsigX_spice
    spice$sigX <- chol2inv(chol(invsigX_spice))
  }
  else{
    sigC = stats::cov(cbind(X,Y))
    sigX = stats::cov(X)
    
    invsigC=r.glasso.std(sigC,lam = lambda_spice)
    
    spice$invsigXcY <- invsigC[1:p,1:p]
    spice$sigXcY <- chol2inv(chol(spice$invsigXcY))
    
    invsigX_spice= r.glasso.std(sigX,lam = lambda_spice)
    spice$invsigX = invsigX_spice
    spice$sigX <- chol2inv(chol(invsigX_spice))
  }
  return(spice)
}