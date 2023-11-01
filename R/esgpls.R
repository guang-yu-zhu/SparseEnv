esgpls <- function(X, Y, u, family='logistic',lambda1=NULL,lambda2=NULL,eps=1e-4,init=NULL,asy=0,verbose=0,method="L-BFGS-B")
{
  #lambda1 <- exp(seq(log(1),log(1e-3),len=5));lambda2 <- exp(seq(log(1),log(1e-3),len=5)); eps2=1e-4;init=NULL;asy=1;verbose=1;method="L-BFGS-B"
  t1 = proc.time()
  X=as.matrix(X)
  Y=as.matrix(Y)
  p = ncol(X)
  r = ncol(Y)
  n = nrow(X)
  if(u==0||u==p) out=egpls(X,Y,u,asy=asy,family=family)
  else{
    if(missing(lambda1)) lambda1 <- exp(seq(log(1),log(1e-3),len=5))
    if(missing(lambda2)) lambda2 <- exp(seq(log(1),log(1e-3),len=5))
    if(missing(init)){
      init_tmp=egpls(X,Y,u,family = family)
      init=init_tmp$Gamma
    }
    GEidx = GE(init)	
    newX = X[, GEidx,drop=FALSE]	
    newinit=init[GEidx,,drop=FALSE]
    
    m1 <- LL_esgpls(newX,Y, u, family=family,weight = rep(1,p-u),
                    lambda=lambda1,eps=eps,init=newinit,verbose=verbose,method=method)    #calculating weight     
    #calculating weight
    A <- m1$G[-(1:u),]
    w_norm <- 1/(rowSums(A^2))
    w_norm[w_norm>=1000]=1000
    #w_norm2 <- 1/(rowSums(Gammahat^2)^2)[(u+1):p]
    if(verbose) {print(m1$Gamma);cat('weights:',w_norm,'\n')}       
    #adaptive lasso step     
    m2 <- LL_esgpls(newX, Y, u, family=family,weight = w_norm, lambda=lambda2,  eps=eps, init=newinit,verbose=verbose,method=method)    #calculating weight     
    out=m2
    out$beta=m2$beta[order(GEidx),,drop=FALSE]
    out$Gamma=m2$Gamma[order(GEidx),,drop=FALSE]
    out$Gamma0=m2$Gamma0[order(GEidx),,drop=FALSE]     
    out$where1=     sort(GEidx[m2$where1])     
    out$where0= setdiff(1:p,GEidx[m2$where1])   
    lambdas = c(m1$lambda,m2$lambda)
    names(lambdas)=c('first stage','second stage')
    out$lambda= lambdas
    out$BIC_seq = NULL
    out$BIC_seq1=m1$BIC_seq
    out$BIC_seq2=m2$BIC_seq
    if(asy){
      sigX=cov(X)*n/(n-1)
      q = out$q
      idx=out$where1
      Gamma=out$Gamma
      beta=out$beta
      alpha=out$alpha
      eta=out$eta
      theta = alpha+X%*%beta;
      idx=out$where1
      idx2=out$where0
      
      tmp = egpls_cov(X,Y,alpha,beta,family = family)
      sigXw = tmp$M
      if(family=='logistic')
        wts = 1/(2+exp(-theta)+exp(theta))
      else if(family=='poisson'){
        wts = c(exp(theta))
      }
      covMatrix0 = 1/mean(wts)*solve(sigXw)
      asyFm = matrix(sqrt(diag(covMatrix0)), p, r);
      Omega=t(Gamma) %*% sigX %*% Gamma;
      if(q>u){
        Gamma_work=Gamma[idx,]
        Gamma0_work = grams(nulbasis(t(Gamma_work)))
        Omega0_work = t(Gamma0_work) %*% sigX[idx, idx] %*% Gamma0_work;
        covMatrix0_work = covMatrix0[idx,idx]
        temp = kronecker(eta,t(Gamma0_work))%*%solve(covMatrix0_work)%*%kronecker(t(eta),Gamma0_work) + kronecker(Omega, solve(Omega0_work)) + kronecker(solve(Omega), Omega0_work)- 2 * kronecker(diag(u), diag(q - u))
        covMatrix = (Gamma_work%*%t(Gamma_work))%*%covMatrix0_work%*%(Gamma_work%*%t(Gamma_work)) + kronecker(t(eta), Gamma0_work) %*% solve(temp) %*% kronecker(eta, t(Gamma0_work))
        asySE_work = matrix(sqrt(diag(covMatrix)), q, r);
        asySE=matrix(0,p,r);
        asySE[idx,]=asySE_work;
        out$covMatrix = covMatrix
        out$asySE = asySE
        #out$asyFm = asyFm
        out$ratio = asyFm / asySE
      } #q>u
      else{
        covMatrix=matrix(0,p,p)  #p by p
        covMatrix[idx,idx]=covMatrix0[idx,idx];
        asySE = matrix(sqrt(diag(covMatrix)), p, r); #p by r
        out$covMatrix = covMatrix
        out$asySE = asySE
        #out$asyFm = asyFm
        out$ratio = asyFm / asySE
      } #q=u
      
    } # asy==1
  }
  
  out$fit_time=as.numeric((proc.time()-t1)[3])
  class(out)='esgpls'
  out
}
