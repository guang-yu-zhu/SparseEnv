spxenvbase <- function(X, Y, u,lambda, weight=NULL, ftol=1e-2, maxiter=1e2,eps2=1e-4,  init=NULL,verbose=verbose) {
  X =as.matrix(X)
  Y =as.matrix(Y)
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)
  if(missing(weight)) weight=rep(1,p-u)
  mX = colMeans(X)
  mY = colMeans(Y)
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))
  
  sigY <- stats::cov(Yc) * (n-1)/n
  sigX <- stats::cov(Xc) * (n-1)/n
  sigXY <- stats::cov(Xc, Yc) * (n-1)/n
  
  tmp = chol(sigX) # t(tmp)%*%tmp =sigX
  invsigX <- chol2inv(tmp) # invsigX=tmp2%*%t(tmp2)
  tmp2 <- backsolve(tmp, diag(p)) # tmp2=inv(tmp)  
  tmp3 <- crossprod(sigXY, tmp2)
  sigYcX = sigY - tcrossprod(tmp3,tmp3) #sigYX %*% invsigX%*% t(sigYX)
  
  tmp = chol(sigY) # t(tmp)%*%tmp =sigY
  invsigY <- chol2inv(tmp) # invsigY=tmp2%*%t(tmp2)
  tmp2 <- backsolve(tmp, diag(r)) # tmp2=inv(tmp)  
  tmp3 <- sigXY %*% tmp2
  U <- tcrossprod(tmp3,tmp3) #sigXY %*% invsigY%*% t(sigXY)
  sigXcY <- sigX - U
  
  betaOLS <- invsigX%*%sigXY
  logDetsigY = logdet(sigY)  
  logDetsigX = logdet(sigX)
  ModelOutput=list()  
  if (u == 0) {
    ModelOutput$mu = mY
    ModelOutput$beta = matrix(0,p, r)
    ModelOutput$Gamma = NULL
    ModelOutput$Gamma0 = diag(p)
    ModelOutput$eta = NULL
    ModelOutput$SigX = sigX
    ModelOutput$Omega = NULL
    ModelOutput$Omega0 = sigX
    ModelOutput$sigYcX = sigYcX 
    ModelOutput$loglik = - n * (r + p) / 2 * (1 + log(2 * pi)) - n / 2 *  (logDetsigX + logDetsigY)
    ModelOutput$paramNum = r + p * (p + 1) / 2 + r * (r + 1) / 2
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$p = p
    ModelOutput$u = u
    ModelOutput$q = 0
    ModelOutput$where1 = NULL
    ModelOutput$where0 = 1:p
    ModelOutput$lambda = lambda
  } else if (u == p) {
    beta = betaOLS;
    ModelOutput$mu = mY - t(beta) %*% mX
	ModelOutput$beta = betaOLS
    ModelOutput$Gamma = diag(p)
    ModelOutput$Gamma0 = NULL
    ModelOutput$eta = beta
    ModelOutput$SigX = sigX;
    ModelOutput$Omega = sigX
    ModelOutput$Omega0 = NULL
    ModelOutput$sigYcX = sigYcX 
    ModelOutput$loglik =  -n * (p + r) * (1 + log(2 * pi)) / 2 - n / 2 * (logDetsigY + logdet(sigXcY));
    ModelOutput$paramNum = r + (p + r) * (p + r + 1) / 2;
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$p = p
    ModelOutput$u = u
    ModelOutput$q = p
    ModelOutput$where1 = 1:p
    ModelOutput$where0 = NULL
    ModelOutput$lambda = lambda
  } 
  else{
    if(missing(init)) {init <- initial_value(Y,X,u)}
    Ginit <- init %*% solve(init[1:u, ]) 
    
    obj1 <- logdet(t(init) %*% sigXcY %*% init) + logdet((t(init) %*% invsigX %*% init))
    widx <- is.finite(weight)      
    Ginit.tmp <- as.matrix(Ginit[(u+1):p, ])
    Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
    obj1 <- obj1 + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
    
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
      while (i < maxiter) {        
        res <- spenvlp(b2=drop(t2), b3=drop(t3), 
                       A1=diag(p-1), A2=sigXcY[p, p]*invC1, A3=invsigX[p, p]*invC2, 
                       lambda=lambda, eps=eps2, maxiter=1e2, 
                       weight=weight[p-u], 
                       a_vec_init=drop(Ginit[p,]))        
        old_Ginit <- Ginit[p, ]
        Ginit[p, ] <- res$a_vec        
        
        
        a <- qr(Ginit)
        Gamma <- qr.Q(a)
        Ginit.tmp <- as.matrix(Ginit[(u+1):p, ])
        Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
        obj5 <- logdet(t(Gamma) %*% sigXcY %*% Gamma) + logdet((t(Gamma) %*% invsigX %*% Gamma)) + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
        if (abs(obj1 - obj5) < ftol * abs(obj1)) {
          break
        }  
        else {
          obj1 <- obj5
          i <- i + 1
        }	
        #if(sum((Ginit[p,]-old_Ginit)^2) < eps) break
        #i <- i + 1		
      }
      if(verbose) cat('The number of iterations:',i,'.\n',sep='')
      a <- qr.Q(qr(Ginit), complete = TRUE)
      Gamma <- a[, 1:u]
      Gamma0 <- a[, p]
    } 
    else {
      #obj1 <- logdet(t(init) %*% sigXcY %*% init) + logdet((t(init) %*% invsigX %*% init))
      #widx <- is.finite(weight)      
      #Ginit.tmp <- as.matrix(Ginit[(u+1):p, ])
      #Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
      #obj1 <- obj1 + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
      
      GUG <- crossprod(Ginit, (sigXcY %*% Ginit))	
      GVG <- crossprod(Ginit, (invsigX %*% Ginit))		
      
      t4 <- crossprod(Ginit[(u+1):p,], Ginit[(u+1):p, ]) + diag(u)
      i <- 1
      while (i < maxiter) {
        #print(i+1)
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
          invt4 <- chol2inv(chol(t4))
          
          res <- spenvlp(b2=drop(t2), b3=drop(t3), A1=invt4, A2=sigXcY[j, j]*invC1, A3=invsigX[j, j]*invC2, lambda=lambda, eps=eps2, maxiter=1e2, weight=weight[j-u], a_vec_init=drop(Ginit[j,]))
          
          
          Ginit[j, ] <- res$a_vec
          g <- as.matrix(Ginit[j, ])
          # print(g)
          t4 <- t4 + tcrossprod(g, g)
          GUGt2 <- g + t2
          GUG <- GUG + tcrossprod(GUGt2, GUGt2) * sigXcY[j, j]
          
          GVGt2 <- g + t3
          GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invsigX[j, j] 
        }
        a <- qr(Ginit)
        Gamma <- qr.Q(a)
        Ginit.tmp <- as.matrix(Ginit[(u+1):p, ])
        Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
        obj5 <- logdet(t(Gamma) %*% sigXcY %*% Gamma) + logdet((t(Gamma) %*% invsigX %*% Gamma)) + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
        if (abs(obj1 - obj5) < ftol * abs(obj1)) {
          break
        }  
        else {
          obj1 <- obj5
          i <- i + 1
        }	
      }
      if(verbose) cat('The number of iterations:',i,'.\n',sep='')
    }
    Gamma <- as.matrix(Gamma)
    if(!is.na(sum(Gamma))) idx <- which(rowSums(abs(Gamma))>0)  
    q = length(idx)
    idx_i = setdiff(1:p, idx);
    if(q>u){
      #---Compute the rest of the parameters based on \Gamma---
      Gamma0 = grams(nulbasis(t(Gamma)))
      Omega = crossprod(Gamma,sigX) %*% Gamma
      invOmega = chol2inv(chol(Omega))
      eta = invOmega %*% t(Gamma) %*% sigXY; 
      Omega0 = crossprod(Gamma0,sigX) %*% Gamma0
      SigX = Gamma %*% Omega %*% t(Gamma) + Gamma0 %*% Omega0 %*% t(Gamma0)
      beta = Gamma %*%  eta
      mu = mY - t(beta) %*% mX      
      # calculate likelihhod
      a = logdet(t(Gamma) %*% sigXcY %*% Gamma);
      b = logdet(t(Gamma) %*% invsigX %*% Gamma);
      l = n * (p + r) * (1 + log(2 * pi)) + n * (a + b + logDetsigX + logDetsigY)
      paramNum = r + u * (r - p + q) + p * (p + 1) / 2 + r * (r + 1) / 2
    }
    else{
      sigX_work=sigX[idx,idx]
      invsigX_work=chol2inv(chol(sigX_work))
      sigXY_work=sigXY[idx,,drop=FALSE]
      sigYcX = sigY - t(sigXY_work) %*% invsigX_work %*% sigXY_work
      Gamma_work = diag(u);
      Gamma = matrix(0,p,u)
      Gamma[idx,] = Gamma_work;        
      Gamma0 = grams(nulbasis(t(Gamma)))
      Omega = t(Gamma) %*% sigX %*% Gamma;
      Omega0 = t(Gamma0) %*% sigX %*% Gamma0
      SigX = Gamma %*% Omega %*% t(Gamma) + Gamma0 %*% Omega0 %*% t(Gamma0)
      beta_work= invsigX_work %*% sigXY_work
      eta=beta_work
      beta=matrix(0,p,r)
      beta[idx,]=beta_work      
      mu = mY - t(beta) %*% mX;
      # calculate BIC and likelihhod
      l = n * (p + r) * (1 + log(2 * pi)) + n * (logdet(sigYcX) + logDetsigX)
      paramNum = r + q * (r - p + q) + p * (p + 1) / 2 + r * (r + 1) / 2
    } #q=u
    ModelOutput$mu = mu;
    ModelOutput$beta = beta
    ModelOutput$Gamma = Gamma
    ModelOutput$Gamma0 = Gamma0
    ModelOutput$eta = eta
    ModelOutput$SigX = SigX
    ModelOutput$Omega = Omega
    ModelOutput$Omega0 = Omega0
    ModelOutput$sigYcX = sigYcX
    ModelOutput$loglik = - 0.5 * l
    ModelOutput$paramNum = paramNum
    ModelOutput$n = n
    ModelOutput$q=q
    ModelOutput$iternum=i
    ModelOutput$lambda=lambda    
    ModelOutput$where1=idx
    ModelOutput$where0=idx_i  
  } # 1<u<r
  ModelOutput
}





