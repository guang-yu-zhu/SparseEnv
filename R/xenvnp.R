xenvnp <- function(X, Y, u,maxiter = 1e2,eps=1e-2,init=NULL,asy=0,verbose=0,spice=null,lambda_spice=0.1) {
  t1 = proc.time()
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)
  mX = colMeans(X)
  mY = colMeans(Y)
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))

  
  if(missing(spice)){
    spice<-spxenv_spice(X,Y,lambda=lambda_spice)
  }
  sigXY <- stats::cov(Xc, Yc) * (n-1)/n
  invsigX=spice$invsigX 
  sigX=spice$sigX
  sigXcY=spice$sigXcY
  invsigXcY=spice$invsigXcY
  tmp = chol(sigX) 
  invsigXL <- backsolve(tmp, diag(p)) 
  betaOLS <- invsigX%*%sigXY
  
  ModelOutput=list()	
  
  if(u == 0){
    ModelOutput$mu = mY;
    ModelOutput$beta = matrix(0,p, r)
    ModelOutput$Gamma = NULL
    ModelOutput$Gamma0 = diag(p);
    ModelOutput$eta = NULL;
    ModelOutput$SigX = sigX
    ModelOutput$Omega = NULL;
    ModelOutput$Omega0 = sigX;	    
    #ModelOutput$loglik = - n * (r + p) / 2 * (1 + log(2 * pi)) - n / 2 *  (logDetsigX + logDetsigY);
    ModelOutput$paramNum = r + p * (p + 1) / 2 + r * (r + 1) / 2;
    ModelOutput$n = n;    
    if(asy==1){
      ModelOutput$covMatrix = NULL;
      ModelOutput$asySE = NULL;
      ModelOutput$ratio = matrix(1,p, r);
    }
  }
  else if (u == p){
    beta = betaOLS
    ModelOutput$mu = mY - t(beta) %*% mX
    ModelOutput$beta = beta;
    ModelOutput$Gamma = diag(p)
    ModelOutput$Gamma0 = NULL
    ModelOutput$eta = beta
    ModelOutput$SigX = sigX;
    ModelOutput$Omega = sigX
    ModelOutput$Omega0 = NULL
    #ModelOutput$loglik =  -n * (p + r) * (1 + log(2 * pi)) / 2 - n / 2 * (logDetsigY + logdet(sigXcY));
    ModelOutput$paramNum = r + (p + r) * (p + r + 1) / 2;
    ModelOutput$n = n
    if(asy==1){
      covMatrix = kronecker(sigYcX, invsigX);
      asyFm = matrix(sqrt(diag(covMatrix)), p,r);  
      ModelOutput$covMatrix = covMatrix
      ModelOutput$asySE = asyFm
      ModelOutput$ratio = matrix(1, p, r)
    }
  }
  else{
    if(missing(init)) init <- initial_value(Y,X,u)
    Ginit <- init %*% solve(init[1:u, ])
    if(u == (p-1)) # now the G is a u-1 by u-1 identity matrix
    {
      U1c2 <- array(0, dim = c(p-1, p-1))
      V1c2 <- array(0, dim = c(p-1, p-1))
      
      U1c2 <- sigXcY[-p, -p] - as.matrix(sigXcY[-p, p]) %*% sigXcY[p, -p] / sigXcY[p, p]
      V1c2 <- invsigX[-p, -p] - as.matrix(invsigX[-p, p]) %*% invsigX[p, -p] / invsigX[p, p]		
      
      t2 <- sigXcY[-p, p] / sigXcY[p, p] # p-1 by 1
      t3 <- invsigX[-p, p] / invsigX[p, p] # p-1 by 1
      invC1 <- chol2inv(chol(U1c2)) # p-1 by p-1
      invC2 <- chol2inv(chol(V1c2)) # p-1 by p-1
      
      fobj <- function(x){
        tmp2 <- x + t2  
        tmp3 <- x + t3
        T2 <- invC1 %*% tmp2  
        T3 <- invC2 %*% tmp3
        -2 * log(1 + sum(x^2)) + log(1 + sigXcY[p, p] * crossprod(tmp2, T2)) + log(1 + invsigX[p, p] * crossprod(tmp3, T3))
      }
      
      gobj <- function(x){
        tmp2 <- x + t2
        tmp3 <- x + t3
        T2 <- invC1 %*% tmp2	
        T3 <- invC2 %*% tmp3
        -4 * x %*% solve(1 + sum(x^2)) + 2 * T2 / as.numeric(1 / sigXcY[p, p] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invsigX[p, p] + crossprod(tmp3, T3))	
      }
      
      i <- 1
      while (i < maxiter) {       
        res <- optim(Ginit[p, ], fobj, gobj, method = "BFGS",control=list(reltol=1e-2))        
        if (abs(fobj(Ginit[p,]) - fobj(res$par)) < ftol * fobj(Ginit[p,])) {
          Ginit[p,] <- res$par
          break
        } else {
          Ginit[p,] <- res$par
          i <- i + 1
        }      
      }
      if(verbose) cat('The number of iterations:',i,'.\n',sep='')
      a <- qr.Q(qr(Ginit), complete = TRUE)
      Gamma <- a[, 1:u]
      Gamma0 <- a[, p]
    } 
	else {    
      GUG <- crossprod(Ginit, (sigXcY %*% Ginit))	
      tmp2 <- crossprod(Ginit,invsigXL)
      GVG <- tcrossprod(tmp2,tmp2)
      #GVG <- crossprod(Ginit, (invsigX %*% Ginit))	
      if(!isSymmetric(GUG)) print('GUG is not symmetric!')
      if(!isSymmetric(GVG)) print('GVG is not symmetric!')
      t4 <- crossprod(Ginit[(u+1):p,], Ginit[(u+1):p, ]) + diag(u)
	  oldGinit=Ginit
      i <- 1
      while (i < maxiter) {      
        for (j in (u+1):p) {
          g <- as.matrix(Ginit[j, ])
          t2 <- crossprod(Ginit[-j, ], as.matrix(sigXcY[-j, j])) / sigXcY[j, j]
          t3 <- crossprod(Ginit[-j, ], as.matrix(invsigX[-j, j])) / invsigX[j, j]
          
          GUGt2 <- g + t2
          GUG <- GUG - tcrossprod(GUGt2, GUGt2) * sigXcY[j, j]
          
          GVGt2 <- g + t3
          GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invsigX[j, j] 
          
          t4 <- t4 - tcrossprod(g, g)
          #save(GUG,GVG,file='GUG.Rdata')          
          invC1 <- chol2inv(sechol(GUG))
          invC2 <- chol2inv(sechol(GVG))
          
          invt4 <- chol2inv(sechol(t4))				
          
          fobj <- function(x) {
            tmp2 <- x + t2
            tmp3 <- x + t3
            T1 <- invt4 %*% x
            T2 <- invC1 %*% tmp2	
            T3 <- invC2 %*% tmp3
            -2 * log(1 + x %*% T1) + log(1 + sigXcY[j, j] * crossprod(tmp2, T2)) + log(1 + invsigX[j, j] * crossprod(tmp3, T3))
          }
          
          gobj <- function(x) {
            tmp2 <- x + t2
            tmp3 <- x + t3
            T1 <- invt4 %*% x
            T2 <- invC1 %*% tmp2	
            T3 <- invC2 %*% tmp3
            -4 	* T1 / as.numeric(1 + x %*% T1) + 2 * T2 / as.numeric(1 / sigXcY[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invsigX[j, j] + crossprod(tmp3, T3))	
          }
          
          res <- optim(Ginit[j,], fobj, gobj, method = "BFGS",control=list(reltol=1e-2))
          Ginit[j, ] <- res$par
          g <- as.matrix(Ginit[j, ])
          t4 <- t4 + tcrossprod(g, g)
          GUGt2 <- g + t2
          GUG <- GUG + tcrossprod(GUGt2, GUGt2) * sigXcY[j, j]
          
          GVGt2 <- g + t3
          GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invsigX[j, j]        
          
        }
        if (norm(Ginit-oldGinit,type='F')<eps){
          break
        }
        else{
          oldGinit = Ginit;
          i = i+1;
        }
      }
      if(verbose) cat('The number of iterations:',i,'.\n',sep='')
      a <- qr(Ginit)
      Gamma <- qr.Q(a)
      Gamma0 <- qr.Q(a, complete = TRUE)[, (u+1):p]      
    }
    
    Omega = crossprod(Gamma,sigX) %*% Gamma
    Omega0 = crossprod(Gamma0,sigX) %*% Gamma0
    SigX = Gamma %*% Omega %*% t(Gamma) + Gamma0 %*% Omega0 %*% t(Gamma0)
    eta = solve(Omega) %*% t(Gamma) %*% sigXY;
    beta = Gamma %*% eta
    mu = mY - t(beta) %*% mX;      
    
    #a = logdet(t(Gamma) %*% sigXcY %*% Gamma);
    #b = logdet(t(Gamma) %*% invsigX %*% Gamma);
    #l = n * (p + r) * (1 + log(2 * pi)) + n * (a + b + logDetsigX + logDetsigY)
    
    ModelOutput$mu = mu
    ModelOutput$beta = beta
    ModelOutput$Gamma = Gamma
    ModelOutput$Gamma0 = Gamma0
    ModelOutput$eta = eta
    ModelOutput$SigX = SigX
    ModelOutput$Omega = Omega
    ModelOutput$Omega0 = Omega0
    #ModelOutput$sigYcX = sigYcX
    #ModelOutput$loglik = - 0.5 * l
    ModelOutput$paramNum = r + u * r + p * (p + 1) / 2 + r * (r + 1) / 2
    ModelOutput$n = n

    
  } 
  class(ModelOutput) <- "xenvnp"
  ModelOutput$fit_time=(proc.time()-t1)[3]
  return(ModelOutput)
}

