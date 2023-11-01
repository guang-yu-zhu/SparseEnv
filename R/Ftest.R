Ftest<-function(X, Y, alpha){
  n = nrow(Y)
  r = ncol(Y)
  p = ncol(X)
  Xc = scale(X,center=TRUE,scale=FALSE)
  Yc = scale(Y,center=TRUE,scale=FALSE)
  pvalue = rep(0, r)
  P = Xc %*% chol2inv(chol(t(Xc) %*% Xc)) %*% t(Xc)
  for (i  in 1:r){
    SSY = crossprod(Yc[, i],Yc[,i])
    Resi = (diag(n) - P) %*% Yc[,i]
    SSE = sum(Resi^2)
    SSM = SSY - SSE
    F = SSM / p / (SSE / (n - 1 - p))
    pvalue[i] = 1 - stats::pf(F, p, n - 1 - p);
  }
  #print(pvalue)
  sort_pv = sort(pvalue)  #*r/(1:r)
  crit_v = alpha / r * (1:r)
  ind_sp = which(sort_pv > crit_v)
  ind_sp
  sort(order(pvalue)[ind_sp])
}

    
