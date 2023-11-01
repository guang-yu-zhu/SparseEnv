init_spxenv <- function(X, Y, u,spice=NULL,lambda_spice=0.1,init_method=1) {
  X <-as.matrix(X)
  Y <-as.matrix(Y)
  p=ncol(X)
  
  if(init_method==1){init=initial_value(Y,X,u)}
  else if(init_method==2) {
	  if(is.null(spice)) {spice=spxenv_spice(X,Y,lambda_spice)}
    init=get_Init(spice$sigXcY,spice$sigX-spice$sigXcY,u)
    }
  #else if(init_method==3)  {init=initial_spls(X,Y,u)}
  #else if(init_method==4) {init = unclass(pls::plsr(Y~X,ncomp=u)$loadings)}
  #else if(init_method==5) {
  #  cv2=glmnet::cv.glmnet(X,Y,nfolds=5)
  #  fit2=glmnet::glmnet(X,Y,lambda = cv2$lambda.min)
  #  where1=which(as.vector(fit2$beta)!=0)
  #  XA=X[,where1]
  #  init = matrix(0,p,u)
  #  tmp = unclass(pls::plsr(Y~XA,ncomp=u)$loadings) #initial_value(Y,XA,u)
  #  init[where1,]=tmp
  #  }	  
  return(init)
}



