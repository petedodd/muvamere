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


##' Simulate responses for single study
##'
##' Simulate responses for single study given covariate data.
##'
##' @title mvn_sample_study
##' @param X covariate data
##' @param betag matrix global regression parameters mean
##' @param sigb matrix global regression parameters SD
##' @param rho interpolant between local and global correlations
##' @param taug global tau mean(s)
##' @param sigt global tau SD(s)
##' @param OmegaG global correlation parameters
##' @param lkj_local local correlation LKJ prior parameter
##' @return matrix of responses
##' @author Pete Dodd
##' @export
##' @import trialr
mvn_sample_study <- function(X, betag, sigb,
                             rho, taug, sigt,
                             OmegaG, lkj_local) {
  ## dimensions
  NV <- ncol(betag)
  Nobs <- nrow(X)
  NP <- ncol(X) # =nrow(betag)
  ## means
  betas <- matrix(rnorm(prod(dim(betag)), mean = c(betag), sd = c(sigb)),
    nrow = nrow(betag), ncol = ncol(betag)
  )
  ## correlations
  taus <- .rtruncnorm0(nrow(OmegaG), mean = taug, sd = sigt)
  OmegaL <- trialr::rlkjcorr(1, K = nrow(OmegaG), eta = lkj_local)
  omega <- (rho * OmegaG + (1 - rho) * OmegaL)
  Sig <- diag(taus) %*% omega %*% diag(taus)
  ## samples
  mvn_simulate(X, betas, Sig)
}

##' Simulate responses for a single, unobserved study using the
##' deviation-penalty (kappa) generative mechanism
##'
##' As \code{mvn_sample_study()}, but for a fitted
##' \code{mvn_infer_mlm_sparse(prior = "rhs_kappa")} model instead of the
##' rho-blend model: the new study's correlations are the global ones plus
##' independent \code{normal(0, kappa)} deviations, resampled (rejection
##' sampling) until positive definite: the same distribution the Stan model
##' itself gives each fitted study, so generation for a genuinely new cohort
##' is consistent with what was fitted. The new study's scales are
##' log-normal, \code{log(tau) ~ N(ltaum, lsig)}, again matching the fitted
##' model. See
##' \code{mvn_generate_AP()}, which dispatches to this or
##' \code{mvn_sample_study()} depending on which model \code{fit} came from.
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
##' generated.
##'
##' @title mvn_simulate_studies
##' @param Xlist list of covariate data matrices
##' @param betag matrix global regression parameters mean
##' @param sigb matrix global regression parameters SD
##' @param rhoA beta prior parameter for interpolant between local and global
##'   correlations
##' @param rhoB beta prior parameter for interpolant between local and global
##'   correlations
##' @param taug global tau mean(s)
##' @param sigt global tau SD(s)
##' @param lkj_local local correlation LKJ prior parameter
##' @param lkj_global global correlation LKJ prior parameter
##' @return a data frame with responses for all studies, with studyno and obsno
##'   columns appended
##' @author Pete Dodd
##' @import trialr
##' @export
mvn_simulate_studies <- function(Xlist,
                                 betag, sigb,
                                 rhoA, rhoB,
                                 taug, sigt,
                                 lkj_local, lkj_global) {
  ## sample globals:
  rho <- rbeta(1, rhoA, rhoB)
  OmegaG <- trialr::rlkjcorr(1, ncol(betag), lkj_global)
  ## loop over studies
  SS <- list()
  for (i in 1:length(Xlist)) {
    SS[[i]] <- mvn_sample_study(
      Xlist[[i]], betag, sigb,
      rho, taug, sigt,
      OmegaG, lkj_local
    )
    SS[[i]] <- as.data.frame(SS[[i]])
    SS[[i]]$studyno <- i
    SS[[i]]$obsno <- 1:nrow(SS[[i]])
  }
  do.call("rbind", SS)
}

## draw from normal(mean, sd) truncated to (0, Inf), matching the rho-blend
## models' fitted prior `tau ~ normal(taum, sigt)` with <lower=0> (a folded
## normal, abs(rnorm()), only agrees with it when mean >> sd). Inverse CDF on
## the log scale, so it stays accurate even when almost all the mass lies
## below zero.
.rtruncnorm0 <- function(n, mean, sd) {
  mean <- rep_len(mean, n)
  sd <- rep_len(sd, n)
  logq <- pnorm(mean / sd, log.p = TRUE) # log P(X > 0)
  v <- log(runif(n)) + logq # log upper-tail prob, uniform on (0, P(X > 0))
  mean - sd * qnorm(v, log.p = TRUE)
}
