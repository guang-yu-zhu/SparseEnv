#' @name Berkeley
#' @title Berkeley guidance data
#' @description
#' Berkeley guidance data (Tuddenham and Snyder, 1953) includes height measurements for 39 boys and 54 girls born in 1928-1929 in Berkeley, CA.
#' @usage
#' data(Berkeley)
#'
#' @format
#' \describe{
#'   \item{X}{X=1 for boy, and X=0 for girl.}
#'   \item{Y}{Height measurements at 31 ages.}
#' }
#'
#' @references
#' Tuddenham, R. D. and M. M. Snyder (1953). Physical growth of California boys and girls from birth
#' to eighteen years. Publications in child development. University of California, Berkeley 1(2), 183–364.
#'
"Berkeley"

#' Pulp and paper property
#'
#' @name fiberpaper
#'
#' @description
#' A data frame with 62 observations on the following 8 variables.
#'
#' @usage
#' data(fiberpaper)
#' @format
#'
#' A data frame with 62 observations on the following 8 variables.
#'
#' \describe{
#'   \item{\code{V1}}{Breaking length.}
#'   \item{\code{V2}}{Elastic modulus.}
#'   \item{\code{V3}}{Stress at failure.}
#'   \item{\code{V4}}{Burst strength.}
#'   \item{\code{V5}}{Arithmetic fiber length.}
#'   \item{\code{V6}}{Long fiber fraction.}
#'   \item{\code{V7}}{Fine fiber fraction.}
#'   \item{\code{V8}}{Zero span tensile.}
#' }
#'
#' @details
#'
#' This data set contains measurements of properties of pulp fibers
#' and the paper made from them.
#'
#' @references
#' Johnson, R.A. and Wichern, D.W. (2007). Applied Multivariate Statistical Analysis, 6th edition.
"fiberpaper"


#' The New Zealand mussels data.
#'
#' @title The New Zealand mussels data
#'
#' @description
#' The New Zealand mussels data.
#'
#' @usage
#' data(mussel)
#'
#' @format
#' The data frame contains the following components:
#'
#' \describe{
#'   \item{X}{X=0 if it was found to have pea crabs living within the shell, X=1 if it was free of pea crabs.}
#'   \item{Y}{Logarithms of shell height, shell width, shell length, shell mass, mussel mass, viscera mass.}
#' }
#'
#'
#' @references
#' Cook, R. Dennis, Bing Li, and Francesca Chiaromonte. "Envelope models for parsimonious and efficient multivariate linear regression." Statist. Sinica 20 (2010): 927-1010.
"mussel"


#' SAT score data.
#'
#' @title SAT score data
#'
#' @description
#' The average SAT score of the fifty states in the U.S. in 1982, as well as six variables that are used to predict the average SAT score.
#'
#' @usage
#' data(SAT)
#'
#' @format
#'
#' The data frame contains the following components:
#'
#' \describe{
#'   \item{sat}{SAT score;}
#'   \item{takers}{percentage of the total eligible students in the state who took the exam;}
#'   \item{income}{the median income of families of the test takers;}
#'   \item{years}{the average number of years that the test takers had formal studies in social sciences, natural sciences, and humanities;}
#'   \item{public}{the percentage of the test takers who attended public secondary schools;}
#'   \item{expend}{the total state expenditure on secondary schools;}
#'   \item{rank}{the median percentile ranking of the test takers within their secondary school classes.}
#' }
#'
#' @source
#' Ramsey, F. and Schafer, D. (2012). The Statistical Sleuth: A Course in Methods of Data Analysis. Boston: Cengage Learning.
#'
#' @examples
#' data(SAT)
#' SAT = as.data.frame(scale(SAT, scale = FALSE))
#' new = model.matrix(sat ~ (takers + income + years + public + expend + rank)^2, SAT)
#' X = as.matrix(new[, c(2:22)])
#' Y = as.matrix(SAT[, 1])
#' m1 = spxenv(X, Y, u = 3)
#'
"SAT"


#' Wheat Protein Data.
#'
#' @title Wheat Protein Data
#'
#' @description
#' The protein content of ground wheat samples.
#'
#' @usage
#' data(wheatprotein)
#'
#' @format
#'
#' A data frame with 50 observations on the following 8 variables:
#'
#' \describe{
#'   \item{\code{V1} to \code{V6}}{Measurements of the reflectance of NIR radiation by the wheat samples at 6 wavelengths in the range 1680-2310 nm. The measurements were made on the log(1/reflectance) scale.}
#'   \item{\code{V7}}{The protein content of each sample (in percent).}
#'   \item{\code{V8}}{Binary indicator, 0 for high protein content and 1 for low protein content. The cutoff point is if the protein content is smaller than 9.75.}
#' }
#'
#' @details
#' The data are the result of an experiment to calibrate a near
#' infrared reflectance (NIR) instrument for measuring the protein
#' content of ground wheat samples. The protein content of each
#' sample (in percent) was measured by the standard Kjeldahl method.
#' In Fearn (1983), the problem is to find a linear combination of
#' the measurements that predicts protein content. The estimated
#' coefficients can then be entered into the instrument allowing the
#' protein content of future samples to be read directly. The first
#' 24 cases were used for calibration, and the last 26 samples were
#' used for prediction.
#'
#' @references
#' Fearn, T. (1983). A misuse of ridge regression in the calibration of a near infrared reflectance instrument.
#'
#' @examples
#' data(wheatprotein)
#' X <- wheatprotein[, 8]
#' Y <- wheatprotein[, 1:6]
#' choose_env(X, Y)
#' m1 = env(X, Y, 1)
#' m1$beta
"wheatprotein"
