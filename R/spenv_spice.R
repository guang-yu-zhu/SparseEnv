spenv_spice<-function(X,Y,lambda=0.1){
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))
  n = nrow(Y)
  sigY <- cov(Yc)*(n-1)/n
  sigX <- cov(Xc)*(n-1)/n 
  sigYX <- cov(Yc, Xc)*(n-1)/n
  tmp = sechol(sigX) # t(tmp)%*%tmp =sigX
  invsigX <- chol2inv(tmp) # invsigX=tmp2%*%t(tmp2)
  betaOLS <- sigYX %*% invsigX
  spice=list()
  invsigY_spice = r.glasso.std(sigY,lam = lambda)
  spice$invsigY = invsigY_spice
  spice$sigY = chol2inv(sechol(invsigY_spice))
  
  
  res <- Yc-Xc%*%t(betaOLS)
  invsigRes_spice <- r.glasso.std(cov(res),lam = lambda)
  spice$sigRes <- chol2inv(chol(invsigRes_spice))
  spice$invsigRes <- invsigRes_spice
  return(spice)
}