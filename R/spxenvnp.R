spxenvnp <- function(X, Y, u, lambda1,lambda2, eps = 1e-2, eps2=1e-4,maxiter = 1e2,init=NULL,spice=NULL,init_method=2, ggg=4,lambda_spice=0.1) {
  #ggg=2;init_method=5;eps = 1e-2; eps2=1e-4;maxiter = 1e2;init=NULL;lambda_spice=0.1
  t1 = proc.time()
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = nrow(Y)
  r = ncol(Y)
  p = ncol(X)
  oldX = X
  
  if(u==0|u==p|lambda1==0) {
    out = xenvnp(X,Y,u)
    out$q=p
    out$where1=1:p
    out$where0=NULL
    class(out) <- "spxenvnp"
    return(out)
  }
  if(missing(lambda2)){
    lambda2=lambda1
  }
  
  if(missing(spice)){
    spice<-spxenv_spice(X,Y,lambda=lambda_spice)
  }
  
  
  if(missing(init)){
    init = init_spxenv(X,Y,u=u,spice=spice,init_method=init_method)
    #subspace(init,G)
  }
  
  
  
  GEidx = GE(init)
  newX = X[, GEidx]	
  newinit = init[GEidx,,drop = FALSE]
  spice$invsigX = spice$invsigX[GEidx,GEidx]
  spice$sigX = spice$sigX[GEidx,GEidx]
  spice$sigXcY = spice$sigXcY[GEidx,GEidx]
  spice$invsigXcY = spice$invsigXcY[GEidx,GEidx]
  spice$sigXY = spice$sigXY[GEidx,,drop=FALSE]
  
  #subspace(newinit,G)
  m1 <- spxenvnpbase(newX, Y, u, eps = eps, eps2=eps2,maxiter = maxiter, lambda=lambda1, weight = rep(1,p-u),init = newinit,spice=spice)
  #m1 <- spxenvbase(newX, Y, u, lambda=lambda, weight = rep(1,p-u),init = newinit,spice=spice)
  #subspace(m1$Gamma[order(GEidx),,drop=FALSE],G)
  #calculating weight
  Gammahat <- m1$Gamma
  w <- Gammahat %*% solve(Gammahat[1:u, ])
  w_norm <- 1/(rowMeans(w^2)^(0.5*ggg))[(u+1):p]
  w_norm
  #print(summary(w_norm))
  #print(w_norm/w_norm2)
  #adaptive lasso step
  #newinit = m1$Gamma # warm start
  m2 <- spxenvnpbase(newX, Y, u, eps = eps,eps2=eps2, maxiter = maxiter, lambda=lambda2, weight = w_norm,init = newinit,spice=spice)
  #subspace(m2$Gamma[order(GEidx),,drop=FALSE],G)
  
  out = m2
  out$Gamma = m2$Gamma[order(GEidx),,drop=FALSE]
  out$beta = m2$beta[order(GEidx),,drop=FALSE]
  out$where1 =  sort(GEidx[m2$where1])
  out$where0 =  setdiff(1:p,out$where1)
  out$iter=c(m1$iter,m2$iter)
  out$diff_seq=NULL
  out$diff_seq1=m1$diff_seq
  out$diff_seq2=m2$diff_seq
  tmp = c(m1$fit_time,m2$fit_time,(proc.time()-t1)[3])
  names(tmp)=c('1st stage','2nd stage','Total')
  out$fit_time = tmp
  class(out)<-'spxenvnp'
  out
}