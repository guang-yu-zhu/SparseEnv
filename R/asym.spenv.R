#' Asymptotic covariance of vec(beta) for a fitted sparse envelope model
#'
#' @usage
#' asym.spenv(object)
#'
#' @param object A fitted spenv object.
#'
#' @return
#' \item{covMatrix}{The asymptotic covariance of vec(beta). An rp by rp matrix. The covariance matrix returned is asymptotic. For the actual standard errors, multiply by 1/n.}
#' \item{asySE}{The asymptotic standard error for elements in beta under the envelope model. An r by p matrix. The standard errors returned are asymptotic, for actual standard errors, multiply by 1/sqrt(n).}
#' \item{ratio}{The asymptotic standard error ratio of the standard multivariate linear regression estimator over the envelope estimator, for each element in beta. An r by p matrix.}
#'
#' @references
#' Su Z, Zhu G, Chen X, Yang Y. Sparse envelope model: efficient estimation and response variable selection in multivariate linear regression. Biometrika. 2016 Sep 1;103(3):579-93.
#'
#' @author Guangyu Zhu <guangyuzhu@uri.edu>
#'
#' @seealso \code{choose_spenv} for choosing the dimension of envelope subspace.
#' @export
asym.spenv<-function(object,...)
{
  u=object$u
  r=object$r
  p=object$p
  q=object$q
  out=list()
  if(u<r&u>0){
    Gamma = object$Gamma
    Gamma0 = object$Gamma0
    eta = object$eta
    Omega = object$Omega
    Omega0 = object$Omega0
    Sigma = object$Sigma
    sigRes = object$sigRes
    sigY = object$sigY
    sigX = object$sigX
    idx = object$where1
    idx_i = object$where0
    idx2 = rep(0,q*p)
    for(k in 1:p){
      idx2[((k-1)*q+1) : (k*q)]= idx+(k-1)*r
    }
    Gamma_work <- Gamma[idx,,drop=FALSE]
    Gamma0_work <-  grams(nulbasis(t(Gamma_work)))

    Sigma1_work= Gamma_work %*% Omega %*% t(Gamma_work)
    Omega0_1 = t(Gamma0_work) %*% Sigma[idx, idx] %*% Gamma0_work; # q-u by q-u
    if(q==r) Omega0_1c2=Omega0_1 # q-u by q-u
    else{Omega0_2 = Sigma[idx_i, idx_i];  # r-q by r-q
    Omega0_12 = t(Gamma0_work) %*% Sigma[idx, idx_i]; # q-u by r-q
    Omega0_1c2 = Omega0_1 - Omega0_12 %*% solve(Omega0_2) %*% t(Omega0_12); # q-u by q-u
    }
    #---compute asymptotic variance and get the ratios---
    asyFm = kronecker(solve(sigX), sigRes);
    asyFm = matrix(sqrt(diag(asyFm)), r, p);
    temp = kronecker(eta %*% sigX %*% t(eta)+Omega, solve(Omega0_1c2))+kronecker(solve(Omega), Omega0_1) - 2 * kronecker(diag(u), diag(q - u))
    covMatrix_work = kronecker(solve(sigX), Sigma1_work) + kronecker(t(eta), Gamma0_work) %*% solve(temp) %*% kronecker(eta, t(Gamma0_work));
    covMatrix = matrix(0,r*p,r*p)
    covMatrix[idx2,idx2] = covMatrix_work
    asySE = matrix(sqrt(diag(covMatrix)), r, p)
    out$covMatrix = covMatrix;
    out$asySE = asySE;
    out$ratio = asyFm / asySE;
  }
  else if (u==0){
    out$covMatrix = NULL;
    out$asySE = NULL;
    out$ratio = matrix(1,r, p);
  }
  else{
    sigX = object$sigX
    sigRes =  object$sigRes
    covMatrix = kronecker(solve(sigX), sigRes)
    asyFm = matrix(sqrt(diag(covMatrix)), r, p)
    out$covMatrix = covMatrix;
    out$asySE = asyFm;
    out$ratio = matrix(1,r, p)    }

  class(out)<-'asym.spenv'
  return(out)
}
