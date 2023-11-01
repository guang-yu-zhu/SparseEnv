LL_esgpls <- function(X, Y, u, family,eps=1e-4,lambda=NULL, weight = NULL,init=NULL,verbose=1,method="L-BFGS-B")
{
  t1 = proc.time()
  X = as.matrix(X)
  Y = as.matrix(Y)
  p = ncol(X)
  if(missing(lambda)) lambda <- exp(seq(log(.1),log(1e-3),len=5)) 
  if(missing(weight)) weight <- rep(1,p-u)
  if(missing(init)){
      init_tmp=egpls_ini(X,Y,u,family=family)
      init=init_tmp$Gamma
    } 
  BIC_seq <- rep(NA, length(lambda))
  loglik_seq <- rep(NA, length(lambda))
  minBIC = Inf
  res=NULL
  idmin = NULL
  for(l in 1:length(lambda)){
    if(verbose) cat('------lambda =',lambda[l],'\n')
    tmp <- esgpls_base(X, Y, u, family=family,lambda=lambda[l], eps=eps,weight=weight,init=init,verbose=verbose,method=method)
    if(verbose) print(tmp$where1)
    BIC <- -2*tmp$loglik +  (tmp$q-u) * u *log(tmp$n)
    BIC_seq[l] <- BIC
    loglik_seq[l]<- -2*tmp$loglik
    if(BIC<minBIC) {minBIC=BIC;res=tmp;idmin=l}
  }
  lambda.min <- lambda[idmin]
  out <- res
  out$lambda=lambda.min  
  if(verbose) {
    tmp=rbind(loglik_seq,BIC_seq) 
    rownames(tmp)=list('loglik','BIC')
    colnames(tmp)=round(lambda,4)
    print(tmp)
    cat('Based on BIC we choose lambda to be',lambda.min,'.\n')
  }
  out$BIC_seq=BIC_seq
  #out$loglik_seq=loglik_seq
  #out$fit_time=(proc.time()-t1)[3]
  out
}