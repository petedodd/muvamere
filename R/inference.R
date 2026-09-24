## this file contains wrappers for stan models

##' MCMC sampling for data from single study
##'
##' This generates a Stan sample for MVN data from a single study
##'
##' @title mvn_infer
##' @param Y multivariate response data - each line sign/symptom values for a
##'   patient
##' @param X multivariable explanatory data - each line a set of predictors for
##'   a patient
##' @param Z optional new covariate data to predict/simulate Y for (NP columns,
##'   same as X; any number of rows). If supplied, the returned stanfit also
##'   carries a generated quantity `Ynew` (NewObs x NV) with one
##'   posterior-predictive draw of Y per MCMC iteration for each row of Z -
##'   extract with \code{mvn_extract_predictions()}. Default NULL fits the
##'   model with no prediction (NewObs=0).
##' @param beta_prior_sd prior for Betas, default=5
##' @param tau_prior_sd = 2.5,
##' @param prior for tau Cauchy scale in correlations, default=2.5
##' @param lkj_prior_scale LKJ prior scale for correlation, default=2
##' @param iter iterations for MCMC, default = 2e3
##' @param cores number of cores to use
##' @param chains number of chains to use
##' @param ... further arguments passed to \code{rstan::sampling()}
##' @return a Stan sample object
##' @author Pete Dodd
##' @import rstan
##' @export
mvn_infer <- function(Y, X, Z = NULL,
                      beta_prior_sd = 5, # prior for Betas, def=5
                      tau_prior_sd = 2.5, # prior for tau,   def=2.5
                      lkj_prior_scale = 2, # prior for cor, def=2
                      iter = 2e3, cores = 4, chains = 4, ...) {
  if (!is.null(Z) && ncol(Z) != ncol(X)) {
    stop(
      "Z must have the same number of columns as X (",
      ncol(X), "), got ", ncol(Z)
    )
  }
  ## prepare data
  sdata <- list(
    Nobs = nrow(Y), # number of observations
    NP = ncol(X), # number of variables
    NV = ncol(Y), # number of variates
    X = X, # covariate data
    Y = Y, # outcomes
    beta_prior_sd = beta_prior_sd, # prior for Betas, def=5
    tau_prior_sd = tau_prior_sd, # prior for tau,   def=2.5
    lkj_prior_scale = lkj_prior_scale, # prior for cor, def=2
    NewObs = if (is.null(Z)) 0 else nrow(Z), # number of new obs to predict for
    Z = if (is.null(Z)) {
      matrix(numeric(0), 0, ncol(X))
    } else {
      Z
    } # covariate data for prediction
  )

  ## sample
  rstan::sampling(stanmodels$mvn_infer1e,
    data = sdata,
    chains = chains,
    cores = cores,
    iter = iter, ...
  )
}



##' Extract posterior-predictive draws for new covariate data from a
##' mvn_infer() fit
##'
##' Pulls the \code{Ynew} generated quantity (populated only if
##' \code{mvn_infer()} was called with a non-NULL \code{Z}) out of a stanfit
##' and reshapes it into a proper 3-D array indexed \[draw, new observation,
##' variate\], rather than leaving the caller to reshape rstan's flattened,
##' row-major summary output by hand (see the row-major reshape note in
##' \code{mvn_extract_hyperparams()}/the package TODO list).
##'
##' @title mvn_extract_predictions
##' @param fit a stanfit object returned by \code{mvn_infer()}, called with
##'   non-NULL Z
##' @return a 3-D array \[ndraws, NewObs, NV\] of posterior-predictive draws of
##'   Ynew
##' @author Pete Dodd
##' @export
mvn_extract_predictions <- function(fit) {
  ## NOTE: calling rstan::extract(fit, pars = "Ynew") directly throws
  ## ("no parameter Ynew") rather than returning NULL when Ynew has zero
  ## rows (NewObs=0)
  all_draws <- rstan::extract(fit)
  if (!("Ynew" %in% names(all_draws)) || prod(dim(all_draws$Ynew)[-1]) == 0) {
    stop("fit has no (non-empty) Ynew: ",
      "was mvn_infer() called with a non-NULL Z?")
  }
  ## rstan::extract() already returns this shaped as [draw, NewObs, NV]
  all_draws$Ynew
}





## The hierarchical Stan model evaluates the likelihood study-by-study and
## needs records in contiguous blocks with study ids 1..Nstudies in order. Sort
## here so callers can pass records in any order; study ids are relabelled
## 1..S in sorted order of their original values.
.sort_by_study <- function(Y, X, study) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  if (nrow(X) != nrow(Y) || length(study) != nrow(Y)) {
    stop("Y, X and study must all have the same number of records")
  }
  o <- order(study)
  list(
    Y = Y[o, , drop = FALSE], X = X[o, , drop = FALSE],
    study = as.integer(factor(study[o], levels = sort(unique(study))))
  )
}


## ---- helpers for the sparse (rhs_kappa) model ----------------

## default prior guess for the number of non-zero correlations: 10% of the
## D_R = choose(NV,2) pairs, rounded up, but always strictly below D_R (needed
## by the Stan model)
.default_p0 <- function(NV) {
  DR <- choose(NV, 2)
  min(ceiling(0.1 * DR), DR / 2)
}

## row-major index of the upper triangle of an NV x NV matrix: this is the
## order in which the Stan model numbers the correlations (loc[i,j], i<j)
.upper_idx <- function(NV) {
  idx <- which(upper.tri(diag(NV)), arr.ind = TRUE)
  idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
}

## pooled within-study residual correlation matrix: OLS of Y on X within each
## study, residuals pooled, then their correlation. Non-finite entries (e.g. a
## constant variate) -> 0.
.pooled_resid_cor <- function(Y, X, study) {
  res <- do.call(rbind, lapply(split(seq_len(nrow(Y)), study), function(ix) {
    stats::lm.fit(X[ix, , drop = FALSE], Y[ix, , drop = FALSE])$residuals
  }))
  r <- suppressWarnings(stats::cor(res))
  r[!is.finite(r)] <- 0
  diag(r) <- 1
  r
}

## Initial values for mvn_inferH_sparse_kappa.stan. The non-centered
## horseshoe needs Omega_global(init) positive definite, so Stan's random
## inits (T ~ 1) fail. Informed start: Omega_global = 0.9 * pooled
## within-study residual correlation (positive definite whenever that is),
## T = 0.3, lam = 1, slab multiplier caux = 1, jittered by N(0, 0.01^2) in z;
## each study's centered correlations rr[i] start AT that Omega_global (so
## every study matrix is positive definite too, whatever kappa is), kappa = 1.
## The log-scale tau hierarchy (tz, ltaum, lsig) takes Stan's random inits.
## A deliberately sparse "alarm" start was needed by the old rho-blend model
## but not by this one (no collapse in 12 adversarial-start fits).
.rhs_kappa_inits <- function(Y, X, study, chains, slab_scale) {
  NV <- ncol(Y)
  DR <- choose(NV, 2)
  r <- .pooled_resid_cor(Y, X, study)
  S <- length(unique(study))
  T0 <- 0.3
  c0 <- slab_scale * sqrt(1) # caux = 1
  lt0 <- c0 / sqrt(c0^2 + T0^2) # lam_tilde at lam = 1
  z_inf <- 0.9 * r[.upper_idx(NV)] / (T0 * lt0)
  lapply(seq_len(chains), function(ch) {
    zg <- z_inf + stats::rnorm(DR, 0, 0.01)
    list(
      ## zg/lam as 1-d arrays: with NV = 2 (DR = 1) a plain length-1 vector
      ## is read by rstan as a scalar and initialisation fails
      zg = array(zg, DR), T = T0, lam = array(1, DR), caux = 1, kappa = 1,
      ## rr as a plain [S, DR] matrix, not a list of vectors (rstan silently
      ## ignores the latter and falls back to random inits)
      rr = matrix(zg * T0 * lt0, S, DR, byrow = TRUE)
    )
  })
}

## Severe-instability check applied after every mvn_infer_mlm_sparse() fit:
## catches badly mixing fits even when the chains agree with each other.
.severely_unstable <- function(diag, div_frac_threshold = 0.2) {
  (is.finite(diag$rhat_max) && diag$rhat_max > 1.5) ||
    (is.finite(diag$ess_min) && diag$ess_min < 5) ||
    (diag$n_divergent / diag$n_draws > div_frac_threshold)
}

## Runs mvn_diagnose() on `fit`, attaches it as attr(fit, "diagnose"), and
## warns if sampling is severely unstable
.attach_stability_check <- function(fit, context) {
  diag <- mvn_diagnose(fit)
  attr(fit, "diagnose") <- diag
  if (.severely_unstable(diag)) {
    warning(
      "severe sampling instability in the returned fit", context,
      ": Rhat max ",
      signif(diag$rhat_max, 4), ", ESS min ", signif(diag$ess_min, 4), ", ",
      diag$n_divergent, "/", diag$n_draws, " divergent transitions. ",
      "Chains can agree with each other while all mixing badly. Treat this ",
      "fit as unreliable; see attr(fit, \"diagnose\") for the full ",
      "mvn_diagnose() output.",
      call. = FALSE
    )
  }
  fit
}


##' MCMC sampling for data from multiple studies, with sparsity-promoting
##' shrinkage on the global correlation structure
##'
##' Fits the hierarchical multivariate regression model used for synthetic
##' population generation: per-study regression coefficients and scales
##' partially pooled across studies, a global correlation matrix
##' \code{Omega_global} with a sparsity-promoting prior, and each study's
##' own correlation matrix deviating from it with an estimated strength
##' \code{kappa} (the "deviation-penalty" model, \code{prior = "rhs_kappa"}).
##'
##' \strong{Model.} \code{Omega_global}'s off-diagonal entries get the
##' \emph{regularized horseshoe} of Piironen and Vehtari (2017,
##' \doi{10.1214/17-EJS1337SI}) in a non-centered parametrization. Each
##' study's correlations are \code{rr[i] ~ normal(Omega_global, kappa)},
##' \code{kappa ~ half-normal(0, kappa_prior_scale)}, with any study
##' correlation matrix that is not positive definite rejected. Regression
##' coefficients are non-centered (\code{Betas[i] = BetaM + BetaS * Bz[i]})
##' and the scales use a non-centered log-normal hierarchy,
##' \code{log(tau[i]) ~ normal(ltaum, lsig)}.
##'
##' \strong{Why this model.} It replaced (2026-09) a "rho-blend" model in
##' which each study's correlation was a blend, via \code{rho}, of
##' \code{Omega_global} and a free per-study correlation matrix. That model
##' was non-identified: any study's correlation could be reproduced by its
##' free local matrix at no prior cost, giving a spurious posterior mode with
##' \code{Omega_global} near zero, exactly the quantity
##' \code{mvn_generate_AP()} needs for new cohorts. The deviation-penalty
##' model recovered \code{Omega_global} 3-4x more accurately and sampled far
##' more stably. Its parametrization (centered deviations, log-normal tau)
##' was then chosen by an ablation: fewer divergences and 1.6-2.7x faster
##' than a non-centered version with identical accuracy, although with weak
##' data (small studies, large \code{kappa}) some divergences can remain.
##'
##' The model needs valid initial values (Stan's random defaults fail, since
##' \code{Omega_global} and the study correlations are built by hand and
##' rejected if not positive definite); unless \code{init} is given, this
##' function starts \code{Omega_global} at 0.9 times the pooled within-study
##' residual correlation and every study's correlations at that value.
##'
##' \strong{Severe-instability check.} Every returned fit is passed through
##' \code{mvn_diagnose()}; the result is stored in
##' \code{attr(fit, "diagnose")}, and a warning is issued if sampling looks
##' severely unstable (\eqn{\hat R>1.5}, ESS\eqn{{}<5}, or more than 20\% of
##' draws divergent).
##'
##' @title mvn_infer_mlm_sparse
##' @param Y multivariate response data - each line sign/symptom values for a
##'   patient
##' @param X multivariable explanatory data - each line a set of predictors for
##'   a patient
##' @param study vector of which study each record belongs to
##' @param betaM_prior_sd prior SD for the global regression coefficient means
##' @param betaS_prior_sd prior SD (half-normal) for the between-study
##'   regression coefficient SDs
##' @param tauM_prior_sd prior SD for the mean of log(tau) across studies
##' @param tauS_prior_sd prior SD (half-normal) for the SD of log(tau) across
##'   studies
##' @param prior the correlation model; only \code{"rhs_kappa"} (the
##'   deviation-penalty model, see Details) is available. Kept as an
##'   argument so that code written for earlier versions still runs.
##' @param p0 prior guess for the number of non-zero correlations among the
##'   \code{choose(ncol(Y),2)} pairs, strictly between 0 and that number.
##'   Sets the scale of the global shrinkage parameter,
##'   \code{p0 / (D - p0) / sqrt(nrow(Y))}, following Piironen and Vehtari's
##'   recipe with \code{1/sqrt(n)} standing in for their \code{sigma/sqrt(n)}
##'   (a heuristic transplant). Default: 10\% of the pairs, rounded up (5 for
##'   10 variates).
##' @param slab_scale scale of the Student-t slab that regularizes large
##'   correlations, default 0.5 (correlations lie in (-1,1))
##' @param slab_df slab degrees of freedom, default 4
##' @param kappa_prior_scale scale of the half-normal prior on \code{kappa},
##'   the shared strength of each study's correlation deviation from
##'   \code{Omega_global}. Default 2 (used throughout validation; not
##'   otherwise tuned).
##' @param init optional Stan \code{init} argument; default \code{NULL} gives
##'   the informed starts described in Details
##' @param iter iterations for MCMC, default
##' @param cores number of cores to use
##' @param chains number of chains to use
##' @param ... further arguments passed to \code{rstan::sampling()}, e.g.
##'   \code{seed}, \code{refresh} or \code{control}. Unless \code{pars} is
##'   given, the per-study \code{Omega_s} and \code{Sigs} are not saved
##'   (they are derived, unused downstream, and double the fit's size); pass
##'   \code{pars = NA} to save every parameter.
##' @return a Stan sample object carrying an attribute \code{"diagnose"}, the
##'   \code{mvn_diagnose()} result for this fit (see Details)
##' @author Pete Dodd
##' @references Piironen J, Vehtari A (2017). Sparsity information and
##'   regularization in the horseshoe and other shrinkage priors. Electronic
##'   Journal of Statistics 11(2):5018-5051. \doi{10.1214/17-EJS1337SI}
##' @export
##' @import rstan
mvn_infer_mlm_sparse <- function(Y, X, study,
                                 betaM_prior_sd = 1,
                                 betaS_prior_sd = 0.5,
                                 tauM_prior_sd = 1,
                                 tauS_prior_sd = 0.5,
                                 prior = "rhs_kappa",
                                 p0 = NULL, slab_scale = 0.5, slab_df = 4,
                                 kappa_prior_scale = 2,
                                 init = NULL,
                                 iter = 2e3, cores = 4, chains = 4, ...) {
  ## explicit check, not match.arg(): that partially matches, so the removed
  ## prior = "rhs" would silently become "rhs_kappa"
  if (!identical(prior, "rhs_kappa")) {
    stop("prior must be \"rhs_kappa\"; the rho-blend priors (\"rhs\", ",
      "\"horseshoe\") and mvn_infer_mlm() were removed in 2026-09")
  }
  ## prepare data (records sorted by study, see .sort_by_study)
  srt <- .sort_by_study(Y, X, study)
  Y <- srt$Y
  X <- srt$X
  study <- srt$study
  if (ncol(Y) < 2) {
    stop("the sparse models need at least 2 variates (columns of Y)")
  }
  DR <- choose(ncol(Y), 2)
  if (is.null(p0)) p0 <- .default_p0(ncol(Y))
  if (!(is.numeric(p0) && length(p0) == 1 && p0 > 0 && p0 < DR)) {
    stop("p0 must be a single number with 0 < p0 < choose(ncol(Y), 2) = ",
      DR)
  }
  if (!(slab_scale > 0 && slab_df > 0)) {
    stop("slab_scale and slab_df must be positive")
  }
  if (!(is.numeric(kappa_prior_scale) && length(kappa_prior_scale) == 1 &&
    kappa_prior_scale > 0)) {
    stop("kappa_prior_scale must be a single positive number")
  }
  shdata <- list(
    Nrecords = nrow(Y), # number of records/patients
    Nstudies = length(unique(study)), # number of distinct studies
    study = study, # which study does each record correspond to? [Nrecords]
    NP = ncol(X), # number of variables
    NV = ncol(Y), # number of variates
    X = X, # covariate data [Nrecords,NP]
    Y = Y, # outcomes [Nrecords,NV]
    betaM_prior_sd = betaM_prior_sd, # prior for Betas
    betaS_prior_sd = betaS_prior_sd, # prior for Betas
    tauM_prior_sd = tauM_prior_sd, # prior for log(tau)
    tauS_prior_sd = tauS_prior_sd, # prior for log(tau)
    p0 = p0,
    slab_scale = slab_scale,
    slab_df = slab_df,
    kappa_prior_scale = kappa_prior_scale
  )
  if (is.null(init)) {
    init <- .rhs_kappa_inits(Y, X, study, chains, slab_scale)
  }
  args <- list(...)
  if (is.null(args$pars)) {
    ## per-study Omega_s/Sigs are derived and unused downstream: not saving
    ## them halves the fit's memory and saves up to ~20% sampling time
    args$pars <- c("Omega_s", "Sigs")
    args$include <- FALSE
  }
  fit <- do.call(rstan::sampling, c(list(stanmodels$mvn_inferH_sparse_kappa,
    data = shdata, chains = chains, cores = cores, iter = iter,
    init = init
  ), args))
  .attach_stability_check(fit, "")
}
