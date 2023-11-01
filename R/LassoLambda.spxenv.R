LassoLambda.spxenv <- function(X, Y, u, ftol=1e-2, maxiter=1e2,eps2=1e-4,lambda=NULL, weight = NULL,init=NULL,verbose=1)
{
  #t1 = proc.time()
  X = as.matrix(X)
  Y = as.matrix(Y)
  p = ncol(X)
  if(missing(init))     init=initial_value(Y,X,u)
  
  if(missing(lambda)) lambda <- exp(seq(log(1),log(1e-5),len=15)) 
  if(missing(weight)) weight <- rep(1,p-u)
  
  model_vec <- rep(NA, length(lambda))
  minBIC = Inf
  res=NULL
  idmin = NULL
  for(l in 1:length(lambda)){
    if(verbose) cat('------lambda =',lambda[l],'\n')
    tmp <- spxenvbase(X, Y, u, lambda=lambda[l], ftol=ftol,	maxiter=maxiter, eps2=eps2,weight=weight,init=init,verbose=verbose)
    init = tmp$Gamma
    BIC <- -2*tmp$loglik + log(tmp$n) * (tmp$q-u) * u
    model_vec[l] <- BIC
    if(BIC<minBIC) {minBIC=BIC;res=tmp;idmin=l}
  }
lambda.min <- lambda[idmin]
out <- res
out$lambda=lambda.min  
if(verbose) {
  print(model_vec)
  cat('Based on BIC we choose lambda to be',lambda.min,'.\n')
}
out$BIC_seq=model_vec
#out$fit_time=(proc.time()-t1)[3]
out
}