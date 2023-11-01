spenvnrbase <- function(X, Y, u, lambda, weight=NULL,eps=1e-2,eps2=1e-6, maxiter=1e2, init=NULL, verbose=0,spice=NULL,df=TRUE) {
  t1 = proc.time()
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)
  if(missing(weight)) weight=rep(1,r-u)
  mX = colMeans(X)
  mY = colMeans(Y)
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))
  sigY <- stats::cov(Yc)*(n-1)/n
  sigX <- stats::cov(Xc)*(n-1)/n 
  sigYX <- stats::cov(Yc, Xc)*(n-1)/n
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
    ModelOutput$q = 0
    ModelOutput$where1 = NULL
    ModelOutput$where0 = 1:r
    ModelOutput$lambda = lambda
  } 
  else if (u == r) {
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
    ModelOutput$q = r
    ModelOutput$where1 = 1:r
    ModelOutput$where0 = NULL
    ModelOutput$lambda=lambda
  } 
  else{
    if(missing(spice)) spice<-spenv_spice(X,Y)
    sigRes=spice$sigRes
    invsigY=spice$invsigY
    
    if(missing(init))  init <- initial_value(X,Y,u) #envnr(X,Y,u)$Gamma
    else init = as.matrix(init)
    Ginit <- init %*% solve(init[1:u, ]) 
    
    
    if (u == (r-1)) {     
      U1c2 <- array(0, dim = c(r-1, r-1))
      V1c2 <- array(0, dim = c(r-1, r-1))      
      U1c2 <- sigRes[-r, -r] - as.matrix(sigRes[-r, r]) %*% sigRes[r, -r] / sigRes[r, r]
      V1c2 <- invsigY[-r, -r] - as.matrix(invsigY[-r, r]) %*% invsigY[r, -r] / invsigY[r, r]    
      
      
      t2 <- sigRes[-r, r] / sigRes[r, r]
      t3 <- invsigY[-r, r] / invsigY[r, r]
      invC1 <- chol2inv(sechol(U1c2))
      invC2 <- chol2inv(sechol(V1c2))
      i <- 1    
      while (i < maxiter) {        
        res <- spenvlp(b2=drop(t2), b3=drop(t3), 
                       A1=diag(r-1), A2=sigRes[r, r]*invC1, A3=invsigY[r, r]*invC2, 
                       lambda=lambda, eps=eps2, maxiter=maxiter, 
                       weight=weight[r-u], 
                       a_vec_init=drop(Ginit[r,]))      
        old_Ginit <- Ginit[r, ]
        Ginit[r, ] <- res$a_vec      
        if(sum((Ginit[r,]-old_Ginit)^2) < eps) break
        i <- i + 1    
      }
    } 
    else {        
      GUG <- crossprod(Ginit, (sigRes %*% Ginit))
      GVG <- crossprod(Ginit, (invsigY %*% Ginit))  
      #save(GVG,file='GVG2.RData')
      t4 <- crossprod(Ginit[(u+1):r,], Ginit[(u+1):r, ]) + diag(u)
      i <- 1
      diff_seq = rep(0,maxiter)
      while (i < maxiter) {
        old_Ginit <- Ginit
        for (j in (u+1):r) {
          g <- as.matrix(Ginit[j, ])
          t2 <- crossprod(Ginit[-j, ], as.matrix(sigRes[-j, j])) / sigRes[j, j]
          t3 <- crossprod(Ginit[-j, ], as.matrix(invsigY[-j, j])) / invsigY[j, j]
          
          GUGt2 <- g + t2
          GUG <- GUG - tcrossprod(GUGt2, GUGt2) * sigRes[j, j]
          
          GVGt2 <- g + t3
          GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invsigY[j, j] 
          #save(GVG,GVGt2,file='GVG.RData')
          
          t4 <- t4 - tcrossprod(as.matrix(Ginit[j, ]), as.matrix(Ginit[j, ]))
          invC1 <- chol2inv(sechol(GUG))
          invC2 <- chol2inv(sechol(GVG))
          invt4 <- chol2inv(sechol(t4))
          
          res <- spenvlp(b2=drop(t2), b3=drop(t3), A1=invt4, A2=sigRes[j, j]*invC1, A3=invsigY[j, j]*invC2, lambda=lambda, eps=eps2, maxiter=maxiter, weight=weight[j-u], a_vec_init=drop(Ginit[j,]))
          
          Ginit[j, ] <- res$a_vec
          g <- as.matrix(Ginit[j, ])
          t4 <- t4 + tcrossprod(g, g)
          GUGt2 <- g + t2
          GUG <- GUG + tcrossprod(GUGt2, GUGt2) * sigRes[j, j]
          
          GVGt2 <- g + t3
          GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invsigY[j, j] 
        }   
        diff = norm(Ginit-old_Ginit,type='F')
        diff_seq[i]=diff
        if(diff<eps) break      
        i <- i + 1  
      }
      #print('diff_seq')
      #print(diff_seq[1:i])
      if(verbose) cat('The number of iterations:',i,'.\n',sep='')      
    }
    Gamma <- as.matrix(qr.Q(qr(Ginit)))
    Gamma0 = grams(nulbasis(t(Gamma)))
    if(!is.na(sum(Gamma))) idx <- which((rowSums(abs(Gamma))>0))
    q=length(idx)
    paramNum = r + u * (p - r + q) + r * (r + 1) / 2  
    eta <- crossprod(Gamma, betaOLS)
    beta <- Gamma %*% eta
    alpha = mY - beta %*% mX  
    
    
    ModelOutput$alpha = alpha;
    ModelOutput$beta = beta;
    ModelOutput$Gamma = Gamma;
    ModelOutput$Gamma0 = Gamma0;
    ModelOutput$eta = eta;
    #ModelOutput$sigRes=sigRes   
    #ModelOutput$sigY=sigY
    #ModelOutput$sigX=sigX

    ModelOutput$paramNum = paramNum;
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$u = u
    ModelOutput$p = p
    ModelOutput$q = q
    ModelOutput$where1 = idx
    ModelOutput$where0 = setdiff(1:r,idx)     
    ModelOutput$lambda=lambda
    ModelOutput$diff_seq=diff_seq[1:i]
    ModelOutput$iter=i
  } 
  #if(df){
  #  Yfit <- matrix(1,n,1) %*% t(ModelOutput$alpha) + X %*% t(ModelOutput$beta)
  #  tmp <- cor(Y[,idx],Yfit[,idx])
  #  ModelOutput$df = sum(diag(tmp))
  #}
  ModelOutput$fit_time=(proc.time()-t1)[3]
  class(ModelOutput)<-'spenvnrbase'
  ModelOutput
}



