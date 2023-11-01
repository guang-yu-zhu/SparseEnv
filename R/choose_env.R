#' Choose u for envelope model.
#'
#' @title choose_env
#'
#' @description
#' Select the dimension of the envelope subspace using Bayesian information
#' criterion, Akaike information criterion, and Likelihood ratio testing.
#'
#' @usage choose_env(X, Y)
#'
#' @param X Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.
#' @param Y Multivariate responses. An n by r matrix, r is the number of responses and n is the number of observations. The responses must be continuous variables.
#' @param dims A vector. The dimensions to be chosen from. The default is 0 to r.
#' @param maxiter Maximum number of iterations. Default value: 100.
#' @param ftol Tolerance parameter for F. Default value: 1e-2.
#' @param verbose Flag for print out model fitting process, logical 0 or 1. Default value: 0.
#'
#'
#' @return
#' \item{result}{Dimensions of the envelope subspace chosen by BIC, AIC, LRT(0.05), and LRT(0.01).}
#' \item{detail}{-2*Loglik, BIC, AIC, degrees of freedom, and p-values for LRT test.}
#'
#' @references
#'
#' The codes are implemented based on the algorithm in Section 4.3 of Cook et al (2010).
#'
#' Cook, R. Dennis, Bing Li, and Francesca Chiaromonte. "Envelope models for parsimonious and efficient multivariate linear regression." Statist. Sinica 20 (2010): 927-1010.
#'
# @author Guangyu Zhu and Zhihua Su
#'
#'
#'
#' @examples
#'
#' data(wheatprotein)
#' X <- wheatprotein[, 8]
#' Y <- wheatprotein[, 1:6]
#' choose_env(X, Y)
#' m1 = env(X, Y, 1)
#' m1$beta
#' @export
choose_env <- function(X, Y, dims=NULL, maxiter=1e2,ftol=1e-2,verbose=0){
  n = nrow(Y)
  r = ncol(Y)
  if(is.null(dims)) dims=0:r

  object=matrix(0,length(dims),5)
  dimnames(object)=list(dims,c('-2*Loglik','BIC','AIC','df','p-value'))


  ModelOutput0 = env(X, Y, r,maxiter=maxiter,ftol=ftol,verbose=verbose)
  object[length(dims),1] = - 2 * ModelOutput0$loglik ;
  object[length(dims),2] = - 2 * ModelOutput0$loglik + log(n) * ModelOutput0$paramNum;
  object[length(dims),3] = - 2 * ModelOutput0$loglik + 2 * ModelOutput0$paramNum;
  object[length(dims),5] = 1;

  for(i in 1:(length(dims)-1)){
    if(verbose) cat('Current dimension: ', dims[i], ".\n",seq='')
    ModelOutput = env(X, Y, dims[i], maxiter=maxiter,ftol=ftol,verbose=verbose)
    object[i,1] = - 2 * ModelOutput$loglik
    object[i,2] = - 2 * ModelOutput$loglik + log(n) * ModelOutput$paramNum
    object[i,3] = - 2 * ModelOutput$loglik + 2 * ModelOutput$paramNum
    chisq = (- 2 * ModelOutput$loglik) - (- 2 * ModelOutput0$loglik)
    df = ModelOutput0$paramNum - ModelOutput$paramNum;
    object[i,4] = df;
    object[i,5] = 1-stats::pchisq(chisq, df);
  }
  if(verbose) cat('\n')
  object2=rep(0,4);
  names(object2)=c("BIC","AIC","LRT(0.05)","LRT(0.01)")
  object2[1]=dims[which.min(object[,2])]
  object2[2]=dims[which.min(object[,3])]
  object2[3]=dims[which(object[,5]>0.05)[1]]
  object2[4]=dims[which(object[,5]>0.01)[1]]
  temp=list(result=object2,detail=object)
  class(temp) <- "choose_env"
  return(temp)
}

