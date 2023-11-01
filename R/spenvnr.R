spenvnr <- function(X, Y, u, lambda1,lambda2, eps = 1e-2, eps2=1e-6,maxiter = 1e2,init=NULL,spice=NULL,init_method = 1,idx1=NULL,idx2=NULL,warm=TRUE) {
  t1 = proc.time()
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = nrow(Y)
  r = ncol(Y)
  #if(missing(lambda)) lambda <- exp(seq(log(10),log(1e-10),len = 15))
  if(missing(idx1)) idx1=length(lambda1)
  if(missing(idx2)) idx2=length(lambda2)
  if(u==0|u==r) {out=envnr(X,Y,u)}
  else{
    if(missing(spice)){ spice<-spenv_spice(X,Y) }
    
    if(missing(init)){
      if(init_method==1) init =  initial_value(X,Y,u)
      else if(init_method==2) init =  initial_value2(X,Y,u) 
      else if(init_method==3) init = get_Init(spice$sigRes,spice$sigY-spice$sigRes,u)
      else if(init_method==4) init = get_Init2(spice$sigRes,spice$sigY-spice$sigRes,u)
    }
    
    
    GEidx = GE(init)
    newY = Y[, GEidx]	
    newinit = init[GEidx,,drop = FALSE]
    spice$invsigY = spice$invsigY[GEidx,GEidx]
    spice$sigRes = spice$sigRes[GEidx,GEidx]
    
    
    if(warm){
      init_stage1=newinit
      for(i in 1:idx1){
        stage1=rep(NA,idx1)
        names(stage1)=lambda1[1:idx1]
        m1 <- spenvnrbase(X, newY, u, lambda=lambda1[i], weight = rep(1,r-u),
                          init = init_stage1,spice=spice, eps = eps, eps2=eps2,
                          maxiter = maxiter)
        init_stage1 <- m1$Gamma
        stage1[i]=m1$q
      }
    }
    else{
      m1 <- spenvnrbase(X, newY, u, lambda=lambda1[idx1], weight = rep(1,r-u),
                        init = newinit, spice=spice, eps = eps, eps2=eps2,
                        maxiter = maxiter)
      stage1=m1$q
    }
    
    #calculating weight
    Gammahat <- m1$Gamma
    w <- Gammahat %*% solve(Gammahat[1:u, ])
    w_norm <- 1/(rowSums(w^2)^2)[(u+1):r]   
    
    #adaptive lasso step
    if(warm){
      init_stage2=newinit
      stage2=rep(NA,idx2)
      names(stage2)=lambda2[1:idx2]
      for(i in 1:idx2){
        m2 <- spenvnrbase(X, newY, u,lambda=lambda2[i], weight = w_norm,
                          init = init_stage2, spice=spice, eps = eps,eps2=eps2,
                          maxiter = maxiter)
        stage2[i] = m2$q
        init_stage2=m2$Gamma 
      }
    }
    else{
      m2 <- spenvnrbase(X, newY, u, lambda=lambda2[idx2],weight = w_norm,init = newinit,spice=spice,eps = eps,eps2=eps2, maxiter = maxiter)
      stage2=m2$q
    }
    
    
    out = m2
    out$alpha = m2$alpha[order(GEidx),,drop=FALSE]
    out$Gamma = m2$Gamma[order(GEidx),,drop=FALSE]
    out$beta = m2$beta[order(GEidx),,drop=FALSE]
    out$where1 =  sort(GEidx[m2$where1])
    out$where0 =  setdiff(1:r,out$where1)
    #out$iter=c(m1$iter,m2$iter)
    out$lambda1=lambda1
    out$lambda2=lambda2
    out$idx1=idx1
    out$idx2=idx2
    out$stage1 = stage1
    out$stage2 = stage2
    out$fit_time = as.numeric((proc.time()-t1)[3])
  }
  class(out)<-'spenvnr'
  out
}