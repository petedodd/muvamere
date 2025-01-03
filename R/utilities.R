## this is for functions not relying on Stan



##' MVN sample for a single study
##'
##' MVN sample for a single study using explanatory variables X, regression coefficient matrix Beta, and covariance matrix Sigma
##'
##' @title mvn_simulate
##' @param X multivariable explanatory data - each line a set of predictors for a patient
##' @param Beta matrix of regression coefficients
##' @param Sigma covariance matrix for response
##' @return a matrix of responses Y, with nrow(Y)==nrow(X) and ncol(Y)==ncol(Sigma)
##' @author Pete Dodd
##' @import MASS
##' @export
mvn_simulate <- function(X, Beta, Sigma) {
  MUMatrix <- X %*% Beta
  t(apply(MUMatrix, 1, function(x) MASS::mvrnorm(1, mu = x, Sigma = Sigma)))
}

