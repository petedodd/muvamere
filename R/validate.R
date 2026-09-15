## this file contains simple truth-vs-estimate comparison helpers
##
## These are deliberately simple (RMSE on matched elements). A proper
## divergence measure (e.g. KL divergence) could be considered TODO


##' Compare a true mean vector/matrix against an estimated or simulated one
##'
##' A simple RMSE-based comparison, e.g. for comparing \code{mvn_infer()}'s
##' fitted \code{Beta} (or a generated population's empirical column means)
##' against known simulation truth.
##'
##' @title mvn_compare_means
##' @param true_mu true mean vector/matrix
##' @param est_mu estimated/simulated mean vector/matrix -- same number of
##'   elements as \code{true_mu} (compared element-by-element in storage order)
##' @return a data frame with columns \code{true}, \code{est}, and \code{rmse}
##'   (the overall root-mean-squared error, repeated on every row for easy
##'   access via \code{result$rmse[1]})
##' @author Pete Dodd
##' @export
mvn_compare_means <- function(true_mu, est_mu) {
  true_v <- as.vector(true_mu)
  est_v <- as.vector(est_mu)
  if (length(true_v) != length(est_v)) {
    stop("true_mu and est_mu must have the same number of elements")
  }
  data.frame(true = true_v, est = est_v, rmse = sqrt(mean((true_v - est_v)^2)))
}


##' Compare a true covariance/correlation matrix against an estimated or simulated one
##'
##' A simple RMSE-based comparison over the upper triangle (including the
##' diagonal) of two same-shaped matrices, e.g. for comparing \code{mvn_infer()}'s
##' fitted \code{Sigma} against known simulation truth.
##'
##' @title mvn_compare_cov
##' @param true_Sigma true covariance/correlation matrix
##' @param est_Sigma estimated/simulated covariance/correlation matrix, same dimensions as \code{true_Sigma}
##' @return a data frame with columns \code{true}, \code{est}, and \code{rmse}
##'   (the overall root-mean-squared error over the upper triangle, repeated
##'   on every row for easy access via \code{result$rmse[1]})
##' @author Pete Dodd
##' @export
mvn_compare_cov <- function(true_Sigma, est_Sigma) {
  if (!identical(dim(true_Sigma), dim(est_Sigma))) {
    stop("true_Sigma and est_Sigma must have the same dimensions")
  }
  ut <- upper.tri(true_Sigma, diag = TRUE)
  true_v <- true_Sigma[ut]
  est_v <- est_Sigma[ut]
  data.frame(true = true_v, est = est_v, rmse = sqrt(mean((true_v - est_v)^2)))
}

