envnr <- function(X, Y, u,maxiter = 100,eps=1e-2,init=NULL,spice=NULL) {  
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
  n = nrow(Y)
  sigY <- cov(Yc)*(n-1)/n
  sigX <- cov(Xc)*(n-1)/n 
  sigYX <- cov(Yc, Xc)*(n-1)/n
  tmp = sechol(sigX) # t(tmp)%*%tmp =sigX
  invsigX <- chol2inv(tmp) # invsigX=tmp2%*%t(tmp2)
  betaOLS <- sigYX %*% invsigX
  
  ModelOutput=list()
  
  if (u == 0) {	
    ModelOutput$alpha = mY;
    ModelOutput$beta = matrix(0,r, p)
    ModelOutput$Gamma = NULL
    ModelOutput$Gamma0 = diag(r);
    ModelOutput$eta = NULL;
    ModelOutput$paramNum = r + u * p + r * (r + 1) / 2;
    ModelOutput$n = n;    
    ModelOutput$r = r
    ModelOutput$u = u
    ModelOutput$p = p
    
    
  } else if (u == r) {	
    ModelOutput$alpha = mY - betaOLS %*% mX;
    ModelOutput$beta = betaOLS;
    ModelOutput$Gamma = diag(r);
    ModelOutput$Gamma0 = NULL;
    ModelOutput$eta = betaOLS;
    ModelOutput$paramNum = r + u * p + r * (r + 1) / 2;
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$u = u
    ModelOutput$p = p
    
  } 
  else{
    if(missing(spice)){
      spice <- spenv_spice(X,Y)
      #res <- Yc-Xc%*%t(betaOLS)   
      #tmp1 <- glasso.cv(res,.1)
      #sigRes <- tmp1$sigma	
      #tmp2 <- glasso.cv(Y,.1)
      #invsigY=tmp2$omega
    }
    sigRes <- spice$sigRes
    invsigY <- spice$invsigY
    
    if(missing(init))     init=initial_value(X,Y,u)
    Ginit <- init %*% solve(init[1:u, ])
    if (u == (r-1)) # now the G is a u-1 by u-1 identity matrix
    {      
      # C1
      U1c2 <- array(0, dim = c(r-1, r-1)) 
      # C2
      V1c2 <- array(0, dim = c(r-1, r-1))
      
      U1c2 <- sigRes[-r, -r] - as.matrix(sigRes[-r, r]) %*% sigRes[r, -r] / sigRes[r, r] 
      V1c2 <- invsigY[-r, -r] - as.matrix(invsigY[-r, r]) %*% invsigY[r, -r] / invsigY[r, r]
      
      invC1 <- chol2inv(sechol(U1c2)) # r-1 by r-1
      invC2 <- chol2inv(sechol(V1c2)) # r-1 by r-1
      t2 <- sigRes[-r, r] / sigRes[r, r] # r-1 by 1
      t3 <- invsigY[-r, r] / invsigY[r, r] # r-1 by 1
      fobj <- function(x) {
        tmp2 <- x + t2  
        tmp3 <- x + t3
        T2 <- invC1 %*% tmp2  
        T3 <- invC2 %*% tmp3
        -2 * log(1 + sum(x^2)) + log(1 + sigRes[r, r] * crossprod(tmp2, T2)) + log(1 + invsigY[r, r] * crossprod(tmp3, T3))
      }	
      gobj <- function(x) {
        tmp2 <- x + t2
        tmp3 <- x + t3
        T2 <- invC1 %*% tmp2	
        T3 <- invC2 %*% tmp3
        -4 * x %*% solve(1 + sum(x^2)) + 2 * T2 / as.numeric(1 / sigRes[r, r] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invsigY[r, r] + crossprod(tmp3, T3))	
      }
      
      i <- 1
      while (i < maxiter) {
        res <- optim(Ginit[r, ,drop=FALSE], fobj, gobj, method = "BFGS",control=list(reltol=1e-2))        
        if (norm(Ginit[r,]- res$par,type='F') < eps) {
          Ginit[r,] <- res$par
          break
        } else {
          Ginit[r,] <- res$par
          i <- i + 1
        }        
      }
      cat('The number of iterations:',i,'.\n',sep='')
      a <- qr.Q(qr(Ginit), complete = TRUE)
      Gamma <- a[, 1:u]
      Gamma0 <- a[, r]      
    } 
    ###########################################
    else {
      GUG <- crossprod(Ginit, (sigRes %*% Ginit))	
      GVG <- crossprod(Ginit, (invsigY %*% Ginit))		
      
      t4 <- crossprod(Ginit[(u+1):r,], Ginit[(u+1):r, ]) + diag(u)
      i <- 1
      oldGinit = Ginit;
      while (i < maxiter) {        
        for (j in (u+1):r) {
          g <- as.matrix(Ginit[j, ])
          # t2: U_{22}^{-1} \G^T \U_{12} 
          t2 <- crossprod(Ginit[-j, ], as.matrix(sigRes[-j, j])) / sigRes[j, j]
          # t3: V_{22}^{-1} \G^T \V_{12} 
          t3 <- crossprod(Ginit[-j, ], as.matrix(invsigY[-j, j])) / invsigY[j, j]
          
          # Calculate C1 see page 4 line 5
          GUGt2 <- g + t2
          GUG <- GUG - tcrossprod(GUGt2, GUGt2) * sigRes[j, j]
          # Calculate C2 = GVGt2
          GVGt2 <- g + t3
          GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invsigY[j, j] 
          
          # t4: I+ A_{-r}^TA_{-r}
          t4 <- t4 - tcrossprod(g, g)
          
          invC1 <- chol2inv(sechol(GUG))
          invC2 <- chol2inv(sechol(GVG))
          invt4 <- chol2inv(sechol(t4))				
          
          fobj <- function(x) {
            tmp2 <- x + t2
            tmp3 <- x + t3
            T1 <- invt4 %*% x
            T2 <- invC1 %*% tmp2	
            T3 <- invC2 %*% tmp3
            -2 * log(1 + x %*% T1) + log(1 + sigRes[j, j] * crossprod(tmp2, T2)) +
              log(1 + invsigY[j, j] * crossprod(tmp3, T3))
          }
          
          gobj <- function(x) {
            tmp2 <- x + t2
            tmp3 <- x + t3
            T1 <- invt4 %*% x
            T2 <- invC1 %*% tmp2	
            T3 <- invC2 %*% tmp3
            -4 	* T1 / as.numeric(1 + x %*% T1) + 2 * T2 / as.numeric(1 / sigRes[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invsigY[j, j] + crossprod(tmp3, T3))	
          }
          
          res <- optim(Ginit[j,], fobj, gobj, method = "BFGS",control=list(reltol=1e-2))         
          Ginit[j, ] <- res$par
          g <- as.matrix(Ginit[j, ])
          # update 
          t4 <- t4 + tcrossprod(g, g)
          GUGt2 <- g + t2
          GUG <- GUG + tcrossprod(GUGt2, GUGt2) * sigRes[j, j]
          
          GVGt2 <- g + t3
          GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invsigY[j, j]				          
        } 
        if (norm(Ginit-oldGinit,type='F')<eps){
          break
        }
        else{
          oldGinit = Ginit;
          i = i+1;
        }
      }
      a <- qr(Ginit)
      Gamma <- qr.Q(a) 
      Gamma0 <- qr.Q(a, complete = TRUE)[, (u+1):r,drop=FALSE]
      cat('The number of iterations:',i,'.\n',sep='')
    }
    #---Compute the rest of the parameters based on \Gamma---
    
    eta = t(Gamma) %*% betaOLS
    beta = Gamma %*% eta
    alpha = mY - beta %*% mX  
    Omega = t(Gamma) %*% sigRes %*% Gamma;
    Omega0 = t(Gamma0) %*% sigY %*% Gamma0
    
    ModelOutput$alpha = alpha;
    ModelOutput$beta = beta;
    ModelOutput$Gamma = Gamma;
    ModelOutput$Gamma0 = Gamma0;
    ModelOutput$eta = eta;
    #ModelOutput$Omega=Omega
    #ModelOutput$Omega0=Omega0
    #ModelOutput$sigRes=sigRes   
    #ModelOutput$sigY=sigY
    #ModelOutput$sigX=sigX
    ModelOutput$paramNum = r + u * p + r * (r + 1) / 2;
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$u = u
    ModelOutput$p = p
  }  
  ModelOutput$fit_time=as.numeric((proc.time()-t1)[3])
  class(ModelOutput)<-'envnr'
  ModelOutput
}

