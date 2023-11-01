#' Choose u and lambda for sparse envelope model.
#'
#' @description
#' Select the dimension of the envelope subspace and the tuning parameter using Bayesian information
#' criterion, Akaike information criterion, and Likelihood ratio testing.
#'
#' @usage choose_spenv(X, Y)
#'
#' @param X Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.
#' @param Y Multivariate responses. An n by r matrix, r is the number of responses and n is the number of observations. The responses must be continuous variables.
#' @param dims A vector. The dimensions to be chosen from. The default is 0 to r.
#' @param maxiter Maximum number of iterations. Default value: 100.
#' @param ftol Tolerance parameter for F. Default value: 1e-2.
#' @param verbose Flag for print out model fitting process, logical 0 or 1. Default value: 0.
#'
#' @return
#' \item{result}{Dimensions of the envelope subspace chosen by BIC, AIC, and LRT, and the corresponding tuning parameter chosen for the dimension.}
#' \item{detail}{-2*Loglik, BIC, AIC, degrees of freedom, and p-values for LRT test.}
#'
#' @references
#'
#' The codes are implemented based on the algorithm in Z. Su, G. Zhu, and X. Chen (2015).
#'
#' Su Z, Zhu G, Chen X, Yang Y. Sparse envelope model: efficient estimation and response variable selection in multivariate linear regression. Biometrika. 2016 Sep 1;103(3):579-93.
#'
# @author Zhihua Su and Guangyu Zhu
#'
#'
# @examples
#'
# data(Berkeley)
# X = Berkeley$X
# Y2 = Berkeley$Y[, c(1, 2, 21, 23)]
# choose_spenv(X, Y2)
# m1 = spenv(X, Y2, 1)
#' @export
choose_spenv <- function(X, Y, dims=NULL,lambda1=NULL,lambda2=NULL,maxiter=1e2,ftol=1e-4,verbose=0){
  n = nrow(Y)
  r = ncol(Y)
  if(is.null(dims)) dims <- 0:r
  if(missing(lambda1)) lambda1 <- exp(seq(log(1),log(1e-5),len = 15))
  if(missing(lambda2)) lambda2 <- exp(seq(log(1),log(1e-5),len = 15))

  dims <- unique(c(dims,r))
  object=matrix(0,length(dims),5)
  dimnames(object)=list(dims,c('-2*Loglik','BIC','AIC','df','p-value'))


  ModelOutput0 = env(X, Y, r,maxiter=maxiter,ftol=ftol,verbose=verbose)
  object[length(dims),1] = - 2 * ModelOutput0$loglik ;
  object[length(dims),2] = - 2 * ModelOutput0$loglik + log(n) * ModelOutput0$paramNum;
  object[length(dims),3] = - 2 * ModelOutput0$loglik + 2 * ModelOutput0$paramNum;
  object[length(dims),5] = 1;

  for(i in 1:(length(dims)-1)){
    if(verbose) cat('=Current dimension: ', dims[i], ".\n",seq='')
    ModelOutput = spenv(X, Y, dims[i],lambda1=lambda1,
                        lambda2=lambda2,maxiter=maxiter,ftol=ftol,verbose=verbose)
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
  class(temp) <- "choose_spenv"
  return(temp)
}

