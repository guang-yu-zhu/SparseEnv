#' Asyptotic covariance of vec(beta) for a fitted envelope model.
#'
#' This function computes the asymptotic covariance of vec(beta) and asymptotic standard error for elements in beta for a fitted envelope model.
#'
#' @param object A fitted env object.
#'
#' @return
#' \item{covMatrix}{The asymptotic covariance of vec(beta). An rp by rp matrix. The covariance matrix returned is asymptotic. For the actual standard errors, multiply by 1/n.}
#' \item{asySE}{The asymptotic standard error for elements in beta under the envelope model. An r by p matrix. The standard errors returned are asymptotic, for actual standard errors, multiply by 1/sqrt(n).}
#' \item{ratio}{The asymptotic standard error ratio of the standard multivariate linear regression estimator over the envelope estimator, for each element in beta. An r by p matrix.}
#'
#' @references
#'
#' Cook, R. Dennis, Bing Li, and Francesca Chiaromonte. "Envelope models for parsimonious and efficient multivariate linear regression." \emph{Statist. Sinica} 20 (2010): 927-1010.
#'
#' @seealso
#'
#' \code{choose_env} for choosing the dimension of the envelope subspace.
#' @export
asym.env<-function(object,...)
{
  u=object$u
  r=object$r
  p=object$p
  out=list()
  if(u<r&u>0){
    Gamma = object$Gamma
    Gamma0 = object$Gamma0
    eta = object$eta
    Omega = object$Omega
    Omega0 = object$Omega0
    sigY = object$sigY
    sigX = object$sigX
    sigRes = object$sigRes
    Sigma = object$Sigma
    P = Gamma%*%t(Gamma)
    Sigma1= P%*%Sigma%*% P
    #---compute asymptotic variance and get the ratios---
    asyFm = kronecker(solve(sigX), sigRes);
    asyFm = matrix(sqrt(diag(asyFm)), r, p);
    temp = kronecker(eta %*% sigX %*% t(eta), solve(Omega0))+
      kronecker(Omega, solve(Omega0)) +
      kronecker(solve(Omega), Omega0) - 2 * kronecker(diag(u), diag(r - u))
    covMatrix = kronecker(solve(sigX), Sigma1) + kronecker(t(eta), Gamma0) %*% solve(temp) %*% kronecker(eta, t(Gamma0));
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
    out$ratio = matrix(1,r, p)
  }

  class(out)<-'asym.env'
  return(out)
}
