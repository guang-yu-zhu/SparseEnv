egpls_cov<- function(X,Y,alpha,beta,family){
  n=nrow(X)
  p = ncol(X)
  theta = c(alpha) + X%*%beta
  
  if(family=='logistic')
  {wts = 1/(2+exp(-theta)+exp(theta))
  wts = c(wts/mean(wts))
  Ys= theta+ (Y- 1/(1+exp(-theta)))/wts
}
  else if(family=='poisson')
  {  wts = c(exp(theta))
  wts = wts/mean(wts)
  Ys = theta+(Y-exp(theta))/wts
}
  	
	
  Exw = t(wts)%*%X/n
  tmpX= (X-matrix(1,n,1)%*%Exw)
  sigXw = t(tmpX)%*%diag(wts)%*%tmpX/n;

  
  Eyw = t(wts)%*%Ys/n;
  tmpy= (Ys-matrix(1,n,1)%*%Eyw)
  sigXyw = t(tmpX)%*%diag(wts)%*%tmpy/n;
  Syw = sum(tmpy^2*wts)/n;
  M = sigXw;
  U = sigXyw%*%t(sigXyw)/Syw
  result=list(M=M,U=U,sigXyw=sigXyw)
}