#' Choose u for the envelope-based PLS.
#'
#' @title choose_xenv
#'
#' @description
#' Select the dimension of the envelope subspace using Bayesian information
#' criterion, Akaike information criterion, and Likelihood ratio testing.
#'
#' @usage choose_xenv(X, Y)
#'
#' @param X Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.
#' @param Y Multivariate responses. An n by r matrix, r is the number of responses and n is the number of observations. The responses must be continuous variables.
#' @param dims A vector. The dimensions to be chosen from. The default is 0 to p.
#' @param maxiter Maximum number of iterations. Default value: 100.
#' @param ftol Tolerance parameter for F. Default value: 1e-2.
#' @param verbose Flag for print out model fitting process, logical 0 or 1. Default value: 0.
#'
#' @return
#' \item{result}{Dimensions of the envelope subspace chosen by BIC, AIC, and LRT.}
#' \item{detail}{-2*Loglik, BIC, AIC, degrees of freedom, and p-values for LRT test.}
#'
#' @references
#' The codes are implemented based on the algorithm in 4.5.1 of Cook et al (2012).
#' Cook, R. Dennis, I. S. Helland, and Zhihua Su. "Envelopes and partial least squares regression." Journal of the Royal Statistical Society Series B 20 (2012): 927-1010.
#' The Grassmann manifold optimization step implements the algorithm in sg_min 2.4.3 by Ross Lippert (http://web.mit.edu/~ripper/www.sgmin.html).
#'
# @author Guangyu Zhu and Zhihua Su
#'
#'
#' @examples
#' data(AIS)
#' choose_xenv(AIS$X, AIS$Y)
#' @export
choose_xenv = function(X, Y, dims=NULL, maxiter=1e2,ftol=1e-4,init_method=1){
  n = nrow(X)
  p = ncol(X)
  if(is.null(range)) dims=0:p



  object=matrix(0,length(dims),5)
  dimnames(object)=list(dims,c('-2*Loglik','BIC','AIC','df','p-value'))


  ModelOutput0 = xenv(X, Y, p,maxiter=maxiter,ftol=ftol)
  object[length(dims),1] = - 2 * ModelOutput0$loglik ;
  object[length(dims),2] = - 2 * ModelOutput0$loglik + log(n) * ModelOutput0$paramNum;
  object[length(dims),3] = - 2 * ModelOutput0$loglik + 2 * ModelOutput0$paramNum;
  object[length(dims),5] = 1;

  for(i in 1:(length(dims)-1)){
    cat('Current dimension: ', dims[i], ".\n",seq='')
    ModelOutput = xenv(X, Y, dims[i], maxiter=maxiter,ftol=ftol,init_method=init_method)

    object[i,1] = - 2 * ModelOutput$loglik
    object[i,2] = - 2 * ModelOutput$loglik + log(n) * ModelOutput$paramNum
    object[i,3] = - 2 * ModelOutput$loglik + 2 * ModelOutput$paramNum
    chisq = (- 2 * ModelOutput$loglik) - (- 2 * ModelOutput0$loglik)
    df = ModelOutput0$paramNum - ModelOutput$paramNum;
    object[i,4] = df;
    object[i,5] = 1-stats::pchisq(chisq, df);
  }
  cat('\n')
  object2=rep(0,4);
  names(object2)=c("BIC","AIC","LRT(0.05)","LRT(0.01)")
  object2[1]=dims[which.min(object[,2])]
  object2[2]=dims[which.min(object[,3])]
  object2[3]=dims[which(object[,5]>0.05)[1]]
  object2[4]=dims[which(object[,5]>0.01)[1]]
  temp=list(result=object2,detail=object)
  class(temp) <- "choose_xenv"
  return(temp)
}




