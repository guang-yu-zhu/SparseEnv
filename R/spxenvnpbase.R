spxenvnpbase <- function(X, Y, u, lambda, weight=NULL,eps=1e-2,eps2=1e-4, maxiter=1e2, init=NULL, verbose=0,spice=NULL) {
  t1 = proc.time()
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)
  mX = colMeans(X)
  mY = colMeans(Y)
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  #sigXY <- stats::cov(Xc, Yc)
  
  if(missing(weight)) weight=rep(1,p-u)
  
  ModelOutput=list()
  
  invsigX=spice$invsigX 
  sigX=spice$sigX
  sigXcY=spice$sigXcY
  invsigXcY=spice$invsigXcY
  sigXY <- spice$sigXY
  
  if(missing(init))  init <- init=init_spxenv(X,Y,u=u,lambda_spice=lambda_spice,init_method=init_method)
  else init = as.matrix(init)
  Ginit <- init %*% solve(init[1:u, ]) 

  if (u == (p-1)) {     
    U1c2 <- array(0, dim = c(p-1, p-1))
    V1c2 <- array(0, dim = c(p-1, p-1))      
    U1c2 <- sigXcY[-p, -p] - as.matrix(sigXcY[-p, p]) %*% sigXcY[p, -p] / sigXcY[p, p]
    V1c2 <- invsigX[-p, -p] - as.matrix(invsigX[-p, p]) %*% invsigX[p, -p] / invsigX[p, p]		
    
    
    t2 <- sigXcY[-p, p] / sigXcY[p, p]
    t3 <- invsigX[-p, p] / invsigX[p, p]
    invC1 <- chol2inv(sechol(U1c2))
    invC2 <- chol2inv(sechol(V1c2))
    i <- 1	
    diff_seq = rep(0,maxiter)
    while (i < maxiter) {        
      res <- spenvlp(b2=drop(t2), b3=drop(t3), 
                     A1=diag(p-1), A2=sigXcY[p, p]*invC1, A3=invsigX[p, p]*invC2, 
                     lambda=lambda, eps=eps2, maxiter=maxiter, 
                     weight=weight[p-u], 
                     a_vec_init=drop(Ginit[p,]))      
      old_Ginit <- Ginit[p, ]
      Ginit[p, ] <- res$a_vec      
      diff = norm(Ginit-old_Ginit,type='F')
      diff_seq[i]=diff	
      if(diff < eps) break
      i <- i + 1		
    }
  } 
  else {        
    GUG <- crossprod(Ginit, (sigXcY %*% Ginit))	
    GVG <- crossprod(Ginit, (invsigX %*% Ginit))    
    t4 <- crossprod(Ginit[(u+1):p,], Ginit[(u+1):p, ]) + diag(u)
    i <- 1
    diff_seq = rep(0,maxiter)
    while (i < maxiter) {
      old_Ginit <- Ginit
      for (j in (u+1):p) {
        g <- as.matrix(Ginit[j, ])
        t2 <- crossprod(Ginit[-j, ], as.matrix(sigXcY[-j, j])) / sigXcY[j, j]
        t3 <- crossprod(Ginit[-j, ], as.matrix(invsigX[-j, j])) / invsigX[j, j]
        
        GUGt2 <- g + t2
        GUG <- GUG - tcrossprod(GUGt2, GUGt2) * sigXcY[j, j]
        
        GVGt2 <- g + t3
        GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invsigX[j, j] 
        
        t4 <- t4 - tcrossprod(as.matrix(Ginit[j, ]), as.matrix(Ginit[j, ]))
        invC1 <- chol2inv(sechol(GUG))
        invC2 <- chol2inv(sechol(GVG))
        invt4 <- chol2inv(sechol(t4))
        
        res <- spenvlp(b2=drop(t2), b3=drop(t3), A1=invt4, A2=sigXcY[j, j]*invC1, A3=invsigX[j, j]*invC2, lambda=lambda, eps=eps2, maxiter=maxiter, weight=weight[j-u], a_vec_init=drop(Ginit[j,]))
        
        Ginit[j, ] <- res$a_vec
        g <- as.matrix(Ginit[j, ])
        t4 <- t4 + tcrossprod(g, g)
        GUGt2 <- g + t2
        GUG <- GUG + tcrossprod(GUGt2, GUGt2) * sigXcY[j, j]
        
        GVGt2 <- g + t3
        GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invsigX[j, j] 
      }   
      diff = norm(Ginit-old_Ginit,type='F')
      diff_seq[i]=diff	  
      if(diff<eps) break      
      i <- i + 1  
      if(verbose) cat('The number of iterations:',i,'.\n',sep='')
    }   
    
  }
  
  Gamma <- as.matrix(qr.Q(qr(Ginit)))
  Gamma0 = grams(nulbasis(t(Gamma)))
  #Omega = crossprod(Gamma,sigX) %*% Gamma
  Omega = stats::cov(Xc%*%Gamma)
  invOmega = chol2inv(chol(Omega))
  eta = invOmega %*% stats::cov(X%*%Gamma,Y); 
  beta = Gamma %*%  eta
  mu = mY - t(beta) %*% mX 
  
  if(!is.na(sum(Gamma))) idx <- which((rowSums(abs(Gamma))>0))
  q=length(idx)
  paramNum = r + u * (p - r + q) + r * (r + 1) / 2	
  
  ModelOutput$mu = mu;
  ModelOutput$beta = beta
  ModelOutput$Gamma = Gamma;
  ModelOutput$Gamma0 = Gamma0;
  ModelOutput$eta = eta;
  ModelOutput$n = n
  ModelOutput$r = r
  ModelOutput$u = u
  ModelOutput$p = p
  ModelOutput$q = q  
  ModelOutput$where1 = idx
  ModelOutput$where0 = setdiff(1:p,idx)	    
  ModelOutput$lambda = lambda
  ModelOutput$diff_seq = diff_seq[1:i]
  ModelOutput$iter = i  
  
  ModelOutput$fit_time=(proc.time()-t1)[3]
  class(ModelOutput)<-'spxenvnpbase'
  ModelOutput 
}



