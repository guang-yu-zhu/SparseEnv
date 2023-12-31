#' Fit the envelope model.
#'
#' This function fits the envelope model to the responses and predictors, using the maximum likelihood estimation. When the dimension of the envelope is between 1 and r-1, we implemented the algorithm in Cook et al. (2016). When the dimension is r, then the envelope model degenerates to the standard multivariate linear regression. When the dimension is 0, it means that X and Y are uncorrelated, and the fitting is different.
#'
#'
#' @param X Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.
#'
#' @param Y Multivariate responses. An n by r matrix, r is the number of responses and n is the number of observations. The responses must be continuous variables.
#'
#' @param u Dimension of the envelope. An integer between 0 and r.
#'
#' @param maxiter Maximum number of iterations. Default value: 100.
#'
#' @param ftol Tolerance parameter for F. Default value: 1e-2.
#'
#' @param init The initial value for the envelope subspace. An r by u matrix. Default value is the one generated by function initial_value.
#'
#' @param verbose Flag for print out model fitting process, logical 0 or 1. Default value: 0.
#'
#' @return
#'
#' \item{alpha}{The estimated intercept in the envelope model. An r by 1 vector.}
#' \item{beta}{The envelope estimator of the regression coefficients. An r by p matrix.}
#' \item{Gamma}{The orthogonal basis of the envelope subspace. An r by u semi-orthogonal matrix.}
#' \item{Gamma0}{The orthogonal basis of the complement of the envelope subspace. An r by r-u semi-orthogonal matrix.}
#' \item{eta}{The coordinates of beta with respect to Gamma. A u by p matrix.}
#' \item{Sigma}{The envelope estimator of the error covariance matrix. An r by r matrix.}
#' \item{Omega}{The coordinates of Sigma with respect to Gamma. A u by u matrix.}
#' \item{Omega0}{The coordinates of Sigma with respect to Gamma0. An r-u by r-u matrix.}
#' \item{loglik}{The maximized log likelihood function. A real number.}
#' \item{paramNum}{The number of parameters in the envelope model. A positive integer.}
#' \item{sigRes}{The sample error covariance matrix. A r by r matrix.}
#' \item{sigY}{The sample response covariance matrix. A r by r matrix.}
#' \item{sigX}{The sample predictor covariance matrix. A p by p matrix.}
#' \item{n}{The number of observations in the data. A positive integer.}
#' \item{r}{The number of responses. A nonnegative integer.}
#' \item{u}{Dimension of the envelope. An integer between 0 and r.}
#' \item{p}{The number of predictors. A positive integer.}
#' \item{fit_time}{The time costs for fitting the envelope model.}
#'
#' @references
#'
#' Cook, R. Dennis, Bing Li, and Francesca Chiaromonte. "Envelope models for parsimonious and efficient multivariate linear regression." \emph{Statist. Sinica} 20 (2010): 927-1010.
#'
#' The codes are implemented based on the algorithm in Cook, R. D., Forzani, L. and Su, Z. (2016), A Note on Fast Envelope Estimation. \emph{Journal of Multivariate Analysis.} 150, 42-54.
#'
#' @seealso
#'
#' \code{choose_env} for choosing the dimension of the envelope subspace.
#' @export
env <- function(X, Y, u,init=NULL,ftol = 1e-2,maxiter = 100,verbose=0) {
  t1 = proc.time()
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)
  mX = colMeans(X)
  mY = colMeans(Y)
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))
  sigY <- stats::cov(Yc)*(n-1)/n
  sigX <- stats::cov(Xc)*(n-1)/n
  sigYX <- stats::cov(Yc, Xc)*(n-1)/n
  tmp = sechol(sigX) # t(tmp)%*%tmp =sigX
  invsigX <- chol2inv(tmp) # invsigX=tmp2%*%t(tmp2)
  tmp2 <- backsolve(tmp, diag(p)) # tmp2=inv(tmp)
  tmp3 <- sigYX %*% tmp2
  U <- tcrossprod(tmp3,tmp3) #sigYX %*% invsigX%*% t(sigYX)
  sigRes <- sigY - U
  betaOLS <- sigYX %*% invsigX
  logDetsigY = logdet(sigY)
  ModelOutput=list()

  if (u == 0) {
    ModelOutput$alpha = mY;
    ModelOutput$beta = matrix(0,r, p)
    ModelOutput$Gamma = NULL
    ModelOutput$Gamma0 = diag(r);
    ModelOutput$eta = NULL;
    ModelOutput$Sigma = sigY
    ModelOutput$Omega = NULL;
    ModelOutput$Omega0 = sigY;
    ModelOutput$loglik = - n * r / 2 * (1 + log(2 * pi)) - n / 2 * logDetsigY;
    ModelOutput$sigRes = sigRes
    ModelOutput$sigY = sigY
    ModelOutput$sigX = sigX
    ModelOutput$paramNum = r + u * p + r * (r + 1) / 2;
    ModelOutput$n = n;
    ModelOutput$r = r
    ModelOutput$u = u
    ModelOutput$p = p
  }
  else if (u == r) {
    ModelOutput$alpha = mY - betaOLS %*% mX;
    ModelOutput$beta = betaOLS;
    ModelOutput$Gamma = diag(r);
    ModelOutput$Gamma0 = NULL;
    ModelOutput$eta = betaOLS;
    ModelOutput$Sigma = sigRes;
    ModelOutput$Omega = sigRes;
    ModelOutput$Omega0 = NULL;
    ModelOutput$loglik = - n * r / 2 * (1 + log(2 * pi)) - n / 2 * logdet(sigRes);
    ModelOutput$paramNum = r + u * p + r * (r + 1) / 2;
    ModelOutput$sigRes = sigRes
    ModelOutput$sigY = sigY
    ModelOutput$sigX = sigX
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$u = u
    ModelOutput$p = p
  }
  else{
    invsigY <- chol2inv(sechol(sigY))
    if(missing(init)) init <- initial_value(X,Y,u)
    Ginit <- init %*% solve(init[1:u, ])
    if (u == (r-1)) # now the G is a u-1 by u-1 identity matrix
    {
      # C1
      U1c2 <- array(0, dim = c(r-1, r-1))
      # C2
      V1c2 <- array(0, dim = c(r-1, r-1))

      U1c2 <- sigRes[-r, -r] - as.matrix(sigRes[-r, r]) %*% sigRes[r, -r] / sigRes[r, r]
      V1c2 <- invsigY[-r, -r] - as.matrix(invsigY[-r, r]) %*% invsigY[r, -r] / invsigY[r, r]

      invC1 <- chol2inv(sechol(U1c2)) # r-1 by r-1
      invC2 <- chol2inv(sechol(V1c2)) # r-1 by r-1
      t2 <- sigRes[-r, r] / sigRes[r, r] # r-1 by 1
      t3 <- invsigY[-r, r] / invsigY[r, r] # r-1 by 1
      fobj <- function(x) {
        tmp2 <- x + t2
        tmp3 <- x + t3
        T2 <- invC1 %*% tmp2
        T3 <- invC2 %*% tmp3
        -2 * log(1 + sum(x^2)) + log(1 + sigRes[r, r] * crossprod(tmp2, T2)) + log(1 + invsigY[r, r] * crossprod(tmp3, T3))
      }
      gobj <- function(x) {
        tmp2 <- x + t2
        tmp3 <- x + t3
        T2 <- invC1 %*% tmp2
        T3 <- invC2 %*% tmp3
        -4 * x %*% solve(1 + sum(x^2)) + 2 * T2 / as.numeric(1 / sigRes[r, r] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invsigY[r, r] + crossprod(tmp3, T3))
      }

      i <- 1
      while (i < maxiter) {
        res <- optim(Ginit[r, ], fobj, gobj, method = "BFGS")
        if (abs(fobj(Ginit[r,]) - fobj(res$par)) < ftol * fobj(Ginit[r,])) {
          Ginit[r,] <- res$par
          break
        } else {
          Ginit[r,] <- res$par
          i <- i + 1
        }
      }
      if(verbose==1) cat('The number of iterations:',i,'.\n',sep='')
      a <- qr.Q(qr(Ginit), complete = TRUE)
      Gamma <- a[, 1:u]
      Gamma0 <- a[, r]
      obj5 = logdet(t(Gamma) %*% sigRes %*% Gamma) + logdet((t(Gamma) %*% invsigY %*% Gamma))
    }
    ###########################################
    else {
      obj1 = logdet(t(init) %*% sigRes %*% init)+logdet(t(init) %*% invsigY %*% init)
      GUG <- crossprod(Ginit, (sigRes %*% Ginit))
      GVG <- crossprod(Ginit, (invsigY %*% Ginit))

      t4 <- crossprod(Ginit[(u+1):r,], Ginit[(u+1):r, ]) + diag(u)
      i <- 1
      time1=0
      while (i < maxiter) {
        for (j in (u+1):r) {
          g <- as.matrix(Ginit[j, ])
          # t2: U_{22}^{-1} \G^T \U_{12}
          t2 <- crossprod(Ginit[-j, ], as.matrix(sigRes[-j, j])) / sigRes[j, j]
          # t3: V_{22}^{-1} \G^T \V_{12}
          t3 <- crossprod(Ginit[-j, ], as.matrix(invsigY[-j, j])) / invsigY[j, j]

          # Calculate C1 see page 4 line 5
          GUGt2 <- g + t2
          GUG <- GUG - tcrossprod(GUGt2, GUGt2) * sigRes[j, j]
          # Calculate C2 = GVGt2
          GVGt2 <- g + t3
          GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invsigY[j, j]

          # t4: I+ A_{-r}^TA_{-r}
          t4 <- t4 - tcrossprod(g, g)
          invC1 <- chol2inv(sechol(GUG))
          invC2 <- chol2inv(sechol(GVG))
          invt4 <- chol2inv(sechol(t4))

          fobj <- function(x) {
            tmp2 <- x + t2
            tmp3 <- x + t3
            T1 <- invt4 %*% x
            T2 <- invC1 %*% tmp2
            T3 <- invC2 %*% tmp3
            -2 * log(1 + x %*% T1) + log(1 + sigRes[j, j] * crossprod(tmp2, T2)) +
              log(1 + invsigY[j, j] * crossprod(tmp3, T3))
          }

          gobj <- function(x) {
            tmp2 <- x + t2
            tmp3 <- x + t3
            T1 <- invt4 %*% x
            T2 <- invC1 %*% tmp2
            T3 <- invC2 %*% tmp3
            -4 	* T1 / as.numeric(1 + x %*% T1) + 2 * T2 / as.numeric(1 / sigRes[j, j] + crossprod(tmp2, T2)) + 2 * T3 / as.numeric(1 / invsigY[j, j] + crossprod(tmp3, T3))
          }

          res <- optim(Ginit[j,], fobj, gobj, method = "BFGS")
          #print(time1)

          Ginit[j, ] <- res$par
          g <- as.matrix(Ginit[j, ])
          # update
          t4 <- t4 + tcrossprod(g, g)
          GUGt2 <- g + t2
          GUG <- GUG + tcrossprod(GUGt2, GUGt2) * sigRes[j, j]

          GVGt2 <- g + t3
          GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invsigY[j, j]
        }
        a <- qr(Ginit)
        Gamma <- qr.Q(a)
        obj5 <- logdet(t(Gamma) %*% sigRes %*% Gamma) + logdet((t(Gamma) %*% invsigY %*% Gamma))
        if (abs(obj1/obj5-1) < ftol) {
          Ginit[j,] <- res$par
          break
        } else {
          obj1 <- obj5
          Ginit[j,] <- res$par
          i <- i + 1
        }
      }
      if(verbose==1) cat('The number of iterations:',i,'.\n',sep='')
      Gamma0 <- qr.Q(a, complete = TRUE)[, (u+1):r]

    }
    #---Compute the rest of the parameters based on \Gamma---
    eta = t(Gamma) %*% betaOLS
    beta = Gamma %*% eta
    alpha = mY - beta %*% mX
    Omega = t(Gamma) %*% sigRes %*% Gamma;
    Omega0 = t(Gamma0) %*% sigY %*% Gamma0
    Sigma1 = Gamma %*% Omega %*% t(Gamma)
    Sigma2 = Gamma0 %*% Omega0 %*% t(Gamma0)
    Sigma = Sigma1 + Sigma2

    ModelOutput$alpha = alpha;
    ModelOutput$beta = beta;
    ModelOutput$Gamma = Gamma;
    ModelOutput$Gamma0 = Gamma0;
    ModelOutput$eta = eta;
    ModelOutput$Sigma = Sigma
    ModelOutput$Omega = Omega
    ModelOutput$Omega0 = Omega0
    ModelOutput$loglik = -0.5* (n * r * (1 + log(2 * pi)) + n * (obj5 + logDetsigY))
    ModelOutput$paramNum = r + u * p + r * (r + 1) / 2;
    ModelOutput$sigRes = sigRes
    ModelOutput$sigY = sigY
    ModelOutput$sigX = sigX
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$u = u
    ModelOutput$p = p
  }
  ModelOutput$fit_time=(proc.time()-t1)[3]
  class(ModelOutput)<-'env'
  ModelOutput
}

