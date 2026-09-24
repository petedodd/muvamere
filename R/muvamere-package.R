#' @description Bayesian hierarchical multivariate meta-regression for
#'   generating synthetic patient populations from pooled individual patient
#'   data: per-study regressions partially pooled across studies, with a
#'   sparsity-shrunk global correlation structure that each study's own
#'   correlations deviate from (see \code{mvn_infer_mlm_sparse()}), and
#'   generators for new cohorts (\code{mvn_generate_AP()}).
#' @keywords internal
#' @useDynLib muvamere, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package
#' version 2.32.6. https://mc-stan.org
#'
"_PACKAGE"
