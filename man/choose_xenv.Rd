% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/choose_xenv.R
\name{choose_xenv}
\alias{choose_xenv}
\title{Choose u for the envelope-based PLS.}
\usage{
choose_xenv(X, Y)
}
\arguments{
\item{X}{Predictors. An n by p matrix, p is the number of predictors. The predictors can be univariate or multivariate, discrete or continuous.}

\item{Y}{Multivariate responses. An n by r matrix, r is the number of responses and n is the number of observations. The responses must be continuous variables.}

\item{dims}{A vector. The dimensions to be chosen from. The default is 0 to p.}

\item{maxiter}{Maximum number of iterations. Default value: 100.}

\item{ftol}{Tolerance parameter for F. Default value: 1e-2.}

\item{verbose}{Flag for print out model fitting process, logical 0 or 1. Default value: 0.}
}
\value{
\item{result}{Dimensions of the envelope subspace chosen by BIC, AIC, and LRT.}
\item{detail}{-2*Loglik, BIC, AIC, degrees of freedom, and p-values for LRT test.}
}
\description{
Select the dimension of the envelope subspace using Bayesian information
criterion, Akaike information criterion, and Likelihood ratio testing.
}
\examples{
data(AIS)
choose_xenv(AIS$X, AIS$Y)
}
\references{
The codes are implemented based on the algorithm in 4.5.1 of Cook et al (2012).
Cook, R. Dennis, I. S. Helland, and Zhihua Su. "Envelopes and partial least squares regression." Journal of the Royal Statistical Society Series B 20 (2012): 927-1010.
The Grassmann manifold optimization step implements the algorithm in sg_min 2.4.3 by Ross Lippert (http://web.mit.edu/~ripper/www.sgmin.html).
}
