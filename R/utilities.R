## this is for functions not relying on Stan


##' MVN sample for a single study
##'
##' MVN sample for a single study using explanatory variables X, regression
##' coefficient matrix Beta, and covariance matrix Sigma
##'
##' @title mvn_simulate
##' @param X multivariable explanatory data - each line a set of predictors for
##'   a patient
##' @param Beta matrix of regression coefficients
##' @param Sigma covariance matrix for response
##' @return a matrix of responses Y, with nrow(Y)==nrow(X) and
##'   ncol(Y)==ncol(Sigma)
##' @author Pete Dodd
##' @import MASS
##' @export
mvn_simulate <- function(X, Beta, Sigma) {
  MUMatrix <- X %*% Beta
  t(apply(MUMatrix, 1, function(x) MASS::mvrnorm(1, mu = x, Sigma = Sigma)))
}


##' Simulate responses for a single, unobserved study using the
##' deviation-penalty (kappa) generative mechanism
##'
##' Draws one new study under the model fitted by
##' \code{mvn_infer_mlm_sparse()}: regression coefficients around
##' \code{betag} with SDs \code{sigb}; correlations equal to the global ones
##' plus independent \code{normal(0, kappa)} deviations, resampled
##' (rejection sampling) until positive definite, i.e. the same distribution
##' the Stan model gives each fitted study, so generation for a genuinely new
##' cohort is consistent with what was fitted; and log-normal scales,
##' \code{log(tau) ~ N(ltaum, lsig)}, again matching the fitted model. See
##' \code{mvn_generate_AP()}, which calls this with plug-in estimates from a
##' fit.
##'
##' @title mvn_sample_study_kappa
##' @param X covariate data
##' @param betag matrix global regression parameters mean
##' @param sigb matrix global regression parameters SD
##' @param kappa deviation strength (correlation-scale SD of a new study's
##'   departure from OmegaG)
##' @param ltaum mean(s) of log(tau) across studies (log-normal tau hierarchy,
##'   as fitted by \code{prior = "rhs_kappa"})
##' @param lsig SD(s) of log(tau) across studies
##' @param OmegaG global correlation parameters
##' @param max_tries maximum positive-definiteness rejection-sampling attempts
##'   before erroring (default 200; only relevant for kappa large enough that
##'   most draws are invalid)
##' @return matrix of responses
##' @author Pete Dodd
##' @export
mvn_sample_study_kappa <- function(X, betag, sigb, kappa, ltaum, lsig,
                                   OmegaG, max_tries = 200) {
  ## dimensions
  NV <- ncol(betag)
  ## means
  betas <- matrix(rnorm(prod(dim(betag)), mean = c(betag), sd = c(sigb)),
    nrow = nrow(betag), ncol = ncol(betag)
  )
  ## correlation: OmegaG + kappa * N(0,1) deviations, rejection-sampled for
  ## positive definiteness (the distribution of each study's rr in
  ## inst/stan/mvn_inferH_sparse_kappa.stan)
  ut <- upper.tri(OmegaG)
  omega <- OmegaG
  for (try in seq_len(max_tries)) {
    dev <- kappa * rnorm(sum(ut))
    omega[ut] <- OmegaG[ut] + dev
    ## symmetrise from the upper triangle
    omega[lower.tri(omega)] <- t(omega)[lower.tri(omega)]
    diag(omega) <- 1
    ev <- eigen(omega, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) > 1e-8) break
    if (try == max_tries) {
      stop(
        "mvn_sample_study_kappa: could not draw a positive-definite ",
        "correlation matrix in ", max_tries, " tries; kappa (", kappa,
        ") may be too large relative to OmegaG's own eigenvalues"
      )
    }
  }
  taus <- exp(rnorm(nrow(OmegaG), mean = ltaum, sd = lsig)) # log-normal
  Sig <- diag(taus) %*% omega %*% diag(taus)
  ## samples
  mvn_simulate(X, betas, Sig)
}

##' Simulate responses for a number of studies
##'
##' Given a list of covariate data inputs, a corresponding set of responses is
##' generated under the deviation-penalty model fitted by
##' \code{mvn_infer_mlm_sparse()}: a global correlation matrix is drawn from
##' an LKJ distribution, then each study gets its own correlations (the global
##' ones plus \code{normal(0, kappa)} deviations), regression coefficients
##' and log-normal scales via \code{mvn_sample_study_kappa()}.
##'
##' @title mvn_simulate_studies
##' @param Xlist list of covariate data matrices
##' @param betag matrix global regression parameters mean
##' @param sigb matrix global regression parameters SD
##' @param kappa deviation strength of each study's correlations from the
##'   global correlation matrix
##' @param ltaum mean(s) of log(tau) across studies
##' @param lsig SD(s) of log(tau) across studies
##' @param lkj_global global correlation LKJ prior parameter
##' @return a data frame with responses for all studies, with studyno and obsno
##'   columns appended; the drawn global correlation matrix is attached as
##'   attribute \code{"OmegaG"}
##' @author Pete Dodd
##' @import trialr
##' @export
mvn_simulate_studies <- function(Xlist,
                                 betag, sigb,
                                 kappa, ltaum, lsig,
                                 lkj_global) {
  ## sample globals:
  OmegaG <- trialr::rlkjcorr(1, ncol(betag), lkj_global)
  ## loop over studies
  SS <- list()
  for (i in seq_along(Xlist)) {
    SS[[i]] <- mvn_sample_study_kappa(
      Xlist[[i]], betag, sigb,
      kappa, ltaum, lsig, OmegaG
    )
    SS[[i]] <- as.data.frame(SS[[i]])
    SS[[i]]$studyno <- i
    SS[[i]]$obsno <- 1:nrow(SS[[i]])
  }
  out <- do.call("rbind", SS)
  attr(out, "OmegaG") <- OmegaG
  out
}
