EnvMU<-function(M,U,u){
  p = nrow(U);
  W = matrix(0,p,u)
  tmp=eigen(U);
  W[,1]=tmp$vectors[,1]
  for (k in 2:u){
    Wk = W[,1:(k-1)]
    Ek = M%*%Wk;
    QEk = diag(p) - Ek%*%solve(t(Ek)%*%Ek)%*%t(Ek)
    tmp = eigen(QEk%*%U%*%QEk);
    W[,k] =  tmp$vectors[,1]
  }
  Gamma = grams(Re(W))
  #Gamma0 = nulbasis(t(W))
  return(Gamma)
}