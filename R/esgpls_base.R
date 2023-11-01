esgpls_base <- function(X, Y, u,lambda,family=family, weight, eps=1e-4,init=NULL,asy=0,verbose=0,method="L-BFGS-B") {
  #t1 = proc.time()
  #weight=rep(1,p-u);family='logistic';maxiter = 1e2;ftol=1e-2;eps2=1e-4;init=NULL;asy=0;verbose=0;method="L-BFGS-B"
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)
  mX = colMeans(X)
  mY = colMeans(Y)
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))
  
  sigX <- stats::cov(Xc) * (n-1)/n
  sigXL = chol(sigX) # sigXL%*%t(sigXL) =sigX =crossprod(sigXL)
  #crossprod(sigXL)-sigX
  invsigXL <-t(base::backsolve(sigXL, diag(p)))  #invsigX=crossprod(invsigXL)
  invsigX=crossprod(invsigXL)
  #invsigX-crossprod(invsigXL)   t(invsigXL)%*%invsigXL
  
  ModelOutput=list()	
  if(missing(init)){
    init_tmp=egpls_ini(X,Y,u,family = family)
    init=init_tmp$Gamma
  }
  # calculate other initial value
  Ginit=init%*%solve(init[1:u,])
  A=Ginit[-(1:u),]
  if(family=='logistic')
  {  suppressWarnings(fit0 <-stats::glm(Y~X%*%Ginit,family =binomial(link = "logit"),maxit=5))}
  else if(family=='poisson')
  {  suppressWarnings(fit0 <-stats::glm(Y~X%*%Ginit,family =poisson,maxit=5))}
  
  alpha = c(fit0$coefficients[1])
  eta = matrix(fit0$coefficients[-1],u,1)
  beta=Ginit%*%eta
  tmp = egpls_cov(X,Y,alpha,beta,family=family)
  sigXw=tmp$M
  sigXyw=tmp$sigXyw
  
  flag=1
  iter=1
  while(flag && iter<=5){
    #cat('--------iter',iter,'\n')
    oldA = A
    par0=c(oldA)
    fobj <-function(par){
      A=matrix(par,p-u,u)
      G<-rbind(diag(rep(1,u)),A)
      theta = alpha + X%*%G%*%solve(t(G)%*%sigXw%*%G)%*%t(G)%*%sigXyw;
      if(family=='logistic')
      {  Cn =  -sum(Y*theta) + sum(log(1+exp(theta)))}  
      else if(family=='poisson')
      {  Cn =  -sum(Y*theta) + sum(exp(theta))}
      
      tmp1=sigXL%*%G
      tmp2=invsigXL%*%G  
      Mn=n/2*(logdet(crossprod(tmp1))+logdet(crossprod(tmp2))-2*logdet(crossprod(G)))
      f = 2*(Cn + Mn)/n + lambda*sum(weight*sqrt(rowSums(A^2)))
      f
    }
    
    gobj<-function(par){
      A=matrix(par,p-u,u)
      G<-rbind(diag(rep(1,u)),A)
      GUG=t(G)%*%sigXw%*%G
      invGUG=chol2inv(sechol(GUG))
      GinvGUG=G%*%invGUG 
      GV=t(G)%*%sigXyw  # u by 1
      eta = invGUG%*%GV
      theta = alpha + X%*%G%*%eta;
      if(family=='logistic')
      {dC = -Y + 1/(1+exp(-theta))}
      else if(family=='poisson')
      {dC = -Y + exp(theta)}
      XC = t(X)%*%dC
      P = sigXw%*%G%*%invGUG%*%t(G)
      M1 = XC%*%t(eta) # q by u
      M2 = sigXyw%*%t(XC)%*%G%*% invGUG
      M3 = -P%*%M1; 
      M4 = -P%*%M2; 
      dCn = M1+M2+M3+M4;
      tmp1=sigXL%*%G
      tmp2=invsigXL%*%G
      dMn=n*(sigX%*%G%*%solve(crossprod(tmp1))+invsigX%*%G%*%solve(crossprod(tmp2))-2*G%*%solve(crossprod(G)) )
      dF=2*(dCn+dMn)/n
      m=sqrt(rowSums(A^2))
      m[m==0] = 1e-4
      D = lambda*diag(weight/m)%*%A
      c(dF[-(1:u),]+D)
    }
    #library(numDeriv)
    #grad(fobj,c(AT))-gobj(c(AT))
    #grad(fobj,c(par0))-gobj(c(par0))
    #------
    res <- optim(par0, fobj, gobj, method = method)   #control=list(factr=1e-4)    #control=list(reltol=1e-8)   
    A=matrix(res$par,p-u,u)
    # use subradient method to devide whether i th row of A is 0
    # G<-rbind(diag(rep(1,u)),A)
    # GUG=t(G)%*%sigXw%*%G
    # invGUG=chol2inv(sechol(GUG))
    # GinvGUG=G%*%invGUG
    # GV=t(G)%*%sigXyw  # u by 1
    # eta = invGUG%*%GV
    # theta = alpha + X%*%G%*%eta;
    # dC = -Y + 1/(1+exp(-theta))
    # XC = t(X)%*%dC
    # P = sigXw%*%G%*%invGUG%*%t(G)
    # M1 = XC%*%t(eta) # q by u
    # M2 = sigXyw%*%t(XC)%*%G%*% invGUG
    # M3 = -P%*%M1;
    # M4 = -P%*%M2;
    # dCn = M1+M2+M3+M4;
    # tmp1=sigXL%*%G
    # tmp2=invsigXL%*%G
    # dMn=n*(sigX%*%G%*%solve(crossprod(tmp1))+invsigX%*%G%*%solve(crossprod(tmp2))-2*G%*%solve(crossprod(G)) )
    # dF=2*(dCn+dMn)/n
    # R=dF[-(1:u),]
    # norm_s=sqrt(rowSums(R^2))/(lambda*weight)
    # norm_s
    # setzero=norm_s<1
    norm_s=sqrt(rowSums(A^2))
    setzero=norm_s<1e-3
    A[setzero,] = 0
    q = p-sum(setzero)
    G = rbind(diag(u),A)
    #------
    # update other parameters
    flag=(sum((oldA-A)^2)>eps)
    iter=iter+1
    a <- qr.Q(qr(G), complete = TRUE)
    Gamma <- a[, 1:u]
    if(family=='logistic')
    {  suppressWarnings(fit1 <-stats::glm(Y~X%*%Gamma,family =binomial(link = "logit"),maxit=5))}
    else if(family=='poisson')
    {  suppressWarnings(fit1 <-stats::glm(Y~X%*%Gamma,family =poisson,maxit=5))}
    
    alpha = c(fit1$coefficients[1])
    eta = matrix(fit1$coefficients[-1],u,1)
    beta= Gamma%*%eta
    tmp = egpls_cov(X,Y,alpha,beta,family);
    sigXw=tmp$M
    sigXyw=tmp$sigXyw
  }
  
  Gamma0= qr.Q(qr(Gamma), complete = TRUE)[, (u+1):p] 
  Omega = crossprod(Gamma,sigX) %*% Gamma
  Omega0 = crossprod(Gamma0,sigX) %*% Gamma0
  SigX = Gamma %*% Omega %*% t(Gamma) + Gamma0 %*% Omega0 %*% t(Gamma0)
  
  theta=alpha+X%*%beta
  if(family=='logistic')
  {  Cn =  sum(Y*theta) - sum(log(1+exp(theta)))}
  else if(family=='poisson')
  {  Cn =  sum(Y*theta) - sum(exp(theta))}
  
  Mn= -n/2 * p * (1 + log(2 * pi))-n/2*logdet(SigX)
  loglik = (Cn + Mn)
  
  ModelOutput$alpha = alpha
  ModelOutput$beta = beta
  ModelOutput$Gamma = Gamma
  ModelOutput$G = G
  ModelOutput$Gamma0 = Gamma0
  ModelOutput$eta = eta
  ModelOutput$SigX = SigX
  ModelOutput$Omega = Omega
  ModelOutput$Omega0 = Omega0
  ModelOutput$loglik = loglik
  ModelOutput$paramNum = r + u * r + p * (p + 1) / 2;
  ModelOutput$n = n
  ModelOutput$q =q
  ModelOutput$where1=setdiff(1:p,c(which(setzero)+u))  
  ModelOutput$where0=c(which(setzero)+u)
  #ModelOutput$family=family 
  #ModelOutput$fit_time=(proc.time()-t1)[3]
  return(ModelOutput)
}

