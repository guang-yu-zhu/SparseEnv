egpls_ini<-function (X,Y,u,family='logistic',maxit=5){
  n=nrow(X)
  p = ncol(X)
  if(family=='logistic')
  {suppressWarnings(fit0 <-stats::glm(Y~X,family =binomial(link = "logit"),maxit=5))}
  else if(family=='poisson')
  {suppressWarnings(fit0 <-stats::glm(Y~X,family =poisson,maxit=5))}
  
  a0 = fit0$coefficients[1]
  b0 = fit0$coefficients[-1]
  
  out=list()
  alpha = a0;
  beta = b0;
  
  for (i in 1:maxit){
    tmp = egpls_cov(X,Y,alpha,beta,family=family)
    M=tmp$M
    U=tmp$U
    #if(method==1) 
    G1 = EnvMU(M,U,u)
    #else if(method==2)    G1 = get_Init(M,U,u)
    #else     G1 = get_Init2(M,U,u)
    if(family=='logistic')
    {suppressWarnings(fit1 <-stats::glm(Y~X%*%G1,family =binomial(link = "logit"),maxit=5))}
    else if(family=='poisson')
    {suppressWarnings(fit1 <-stats::glm(Y~X%*%G1,family =poisson,maxit=5))}
    
    alpha = fit1$coefficients[1]
    beta = G1%*%fit1$coefficients[-1];
  }
  G <- G1 %*% solve(G1[1:u, ])
  out$Gamma=G1
  out$G=G
  out
}

