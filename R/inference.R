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
##' @param ...
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
    stop("fit has no (non-empty) Ynew: was mvn_infer() called with a non-NULL Z?")
  }
  ## rstan::extract() already returns this correctly shaped as [draw, NewObs, NV]
  all_draws$Ynew
}





## The hierarchical Stan models evaluate the likelihood study-by-study and need
## records in contiguous blocks with study ids 1..Nstudies in order. Sort here so
## callers can pass records in any order; study ids are relabelled 1..S in
## sorted order of their original values.
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


##' MCMC sampling for data from multiple studies
##'
##' This generates a Stan sample for MVN data from a multiple studies
##'
##' @title mvn_infer_mlm
##' @param Y multivariate response data - each line sign/symptom values for a
##'   patient
##' @param X multivariable explanatory data - each line a set of predictors for
##'   a patient
##' @param study vector of which study each record belongs to
##' @param betaM_prior_sd prior for Betas
##' @param betaS_prior_sd prior for Betas
##' @param tauM_prior_sd prior for tau
##' @param tauS_prior_sd prior for tau
##' @param lkj_local_prior_scale prior for local cor
##' @param lkj_global_prior_scale prior for global cor
##' @param rhoA beta parameter for rho
##' @param rhoB beta parameter for rho
##' @param iter iterations for MCMC, default
##' @param cores number of cores to use
##' @param chains number of chains to use
##' @param ...
##' @return a Stan sample object
##' @author Pete Dodd
##' @export
##' @import rstan
mvn_infer_mlm <- function(Y, X, study,
                          betaM_prior_sd = 1, # prior for Betas
                          betaS_prior_sd = 0.5, # prior for Betas
                          tauM_prior_sd = 1, # prior for tau
                          tauS_prior_sd = 0.5, # prior for tau
                          lkj_local_prior_scale = 3, # prior for local cor
                          lkj_global_prior_scale = 2, # prior for global cor
                          rhoA = 2, # beta parameter for rho
                          rhoB = 2, # beta parameter for rho
                          iter = 2e3, cores = 4, chains = 4, ...) {
  ## prepare data (records sorted by study, see .sort_by_study)
  srt <- .sort_by_study(Y, X, study)
  Y <- srt$Y
  X <- srt$X
  study <- srt$study
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
    tauM_prior_sd = tauM_prior_sd, # prior for tau
    tauS_prior_sd = tauS_prior_sd, # prior for tau
    lkj_local_prior_scale = lkj_local_prior_scale, # prior for local cor
    lkj_global_prior_scale = lkj_global_prior_scale, # prior for global cor
    rhoA = rhoA, # beta parameter for rho
    rhoB = rhoB # beta parameter for rho
  )

  ## sample
  fit <- rstan::sampling(stanmodels$mvn_inferH,
    data = shdata,
    chains = chains,
    cores = cores,
    iter = iter, ...
  )
  .attach_stability_check(fit, "")
}



## ---- helpers for the regularized-horseshoe sparse model ----------------

## default prior guess for the number of non-zero correlations: 10% of the D_R = choose(NV,2)
## pairs, rounded up, but always strictly below D_R (needed by the Stan model)
.default_p0 <- function(NV) {
  DR <- choose(NV, 2)
  min(ceiling(0.1 * DR), DR / 2)
}

## row-major index of the upper triangle of an NV x NV matrix: this is the order in which the
## Stan model numbers the correlations (loc[i,j], i<j)
.upper_idx <- function(NV) {
  idx <- which(upper.tri(diag(NV)), arr.ind = TRUE)
  idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
}

## pooled within-study residual correlation matrix: OLS of Y on X within each study,
## residuals pooled, then their correlation. Non-finite entries (e.g. a constant variate) -> 0.
.pooled_resid_cor <- function(Y, X, study) {
  res <- do.call(rbind, lapply(split(seq_len(nrow(Y)), study), function(ix) {
    stats::lm.fit(X[ix, , drop = FALSE], Y[ix, , drop = FALSE])$residuals
  }))
  r <- suppressWarnings(stats::cor(res))
  r[!is.finite(r)] <- 0
  diag(r) <- 1
  r
}

## Initial values for mvn_inferH_sparse_rhs.stan. The non-centered model needs
## Omega_global(init) positive definite, so Stan's default random inits (T ~ 1) fail.
##   informed : Omega_global = 0.9 * pooled residual correlation (positive definite whenever r is,
##              unit diagonal), T = 0.3, lam = 1, slab multiplier caux = 1, jittered by N(0, 0.01^2) in z
##   sparse   : T = 0.01, z ~ N(0,1) - used only as the alarm start (see check_starts): on a dense
##              truth chains started here can fall into a spurious all-zero mode, which is exactly
##              what disagreement with the informed chains reveals.
## Returns a list of length `chains`; the last chain is the alarm chain if `alarm` is TRUE.
.rhs_inits <- function(Y, X, study, chains, slab_scale, alarm) {
  NV <- ncol(Y)
  DR <- choose(NV, 2)
  r <- .pooled_resid_cor(Y, X, study)
  T0 <- 0.3
  c0 <- slab_scale * sqrt(1) # caux = 1
  lt0 <- c0 / sqrt(c0^2 + T0^2) # lam_tilde at lam = 1
  z_inf <- 0.9 * r[.upper_idx(NV)] / (T0 * lt0)
  lapply(seq_len(chains), function(ch) {
    if (alarm && ch == chains) {
      list(zg = stats::rnorm(DR), T = 0.01, lam = rep(1, DR), caux = 1)
    } else {
      list(zg = z_inf + stats::rnorm(DR, 0, 0.01), T = T0, lam = rep(1, DR), caux = 1)
    }
  })
}

## Agreement between the alarm chain and the other chains: the largest absolute difference in the
## posterior mean of any off-diagonal Omega_global entry (correlation scale). A collapsed alarm chain
## in the dense case differed by ~0.33; agreeing chains by <= 0.005.
.start_gap <- function(fit, alarm_chain) {
  NV <- fit@par_dims$Omega_global[1]
  A <- rstan::extract(fit, pars = "Omega_global", permuted = FALSE) # [iter, chain, NV*NV]
  nch <- dim(A)[2]
  ut <- which(upper.tri(diag(NV))) # column-major positions, matching rstan's flattening
  m <- vapply(seq_len(nch), function(ch) colMeans(A[, ch, ])[ut], numeric(length(ut)))
  others <- rowMeans(m[, -alarm_chain, drop = FALSE])
  max(abs(m[, alarm_chain] - others))
}

## Initial values for mvn_inferH_sparse_kappa.stan (the deviation-penalty
## model): informed Omega_global start as in .rhs_inits, each study's centered
## correlations rr[i] started AT that Omega_global (so every Omega_s[i] is
## positive definite at the start, whatever kappa is), and kappa = 1. The
## log-scale tau hierarchy (tz, ltaum, lsig) takes Stan's random inits.
## Unlike .rhs_inits there is no "alarm" variant as this was found
## empirically unnecessary.
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
      zg = zg, T = T0, lam = rep(1, DR), caux = 1, kappa = 1,
      ## rr as a plain [S, DR] matrix, not a list of vectors (rstan silently
      ## ignores the latter and falls back to random inits)
      rr = matrix(zg * T0 * lt0, S, DR, byrow = TRUE)
    )
  })
}

## Severe-instability check applied uniformly after every mvn_infer_mlm()/mvn_infer_mlm_sparse()
## fit (any prior), regardless of any start-agreement check. For case where chains agree but mix badly.
.severely_unstable <- function(diag, div_frac_threshold = 0.2) {
  (is.finite(diag$rhat_max) && diag$rhat_max > 1.5) ||
    (is.finite(diag$ess_min) && diag$ess_min < 5) ||
    (diag$n_divergent / diag$n_draws > div_frac_threshold)
}

## Runs mvn_diagnose() on `fit`, attaches it as attr(fit, "diagnose"), and warns
.attach_stability_check <- function(fit, context) {
  diag <- mvn_diagnose(fit)
  attr(fit, "diagnose") <- diag
  if (.severely_unstable(diag)) {
    warning(
      "severe sampling instability in the returned fit", context, ": Rhat max ",
      signif(diag$rhat_max, 4), ", ESS min ", signif(diag$ess_min, 4), ", ",
      diag$n_divergent, "/", diag$n_draws, " divergent transitions. ",
      "This is independent of, and not caught by, any start-agreement check: ",
      "two chains can agree with each other while both mixing badly. Treat this fit as ",
      "unreliable; see attr(fit, \"diagnose\") for the full mvn_diagnose() output.",
      call. = FALSE
    )
  }
  fit
}


##' MCMC sampling for data from multiple studies, with sparsity-promoting
##' shrinkage on the global correlation structure
##'
##' As \code{mvn_infer_mlm()}, but the global correlation
##' \code{Omega_global} gets a shrinkage prior on its off-diagonal
##' entries (promoting sparsity as the number of variates grows) instead of
##' a plain LKJ prior. Each study's local correlation \code{Omega_local[i]}
##' uses a plain LKJ prior.
##'
##' Three priors are available. \code{prior = "rhs_kappa"} (the
##' \strong{default}, see below) is a deviation-penalty model that replaces
##' both \code{rho} and \code{Omega_local[i]}'s free per-study prior with a
##' single estimated deviation strength \code{kappa}: each study's
##' correlations are \code{rr[i] ~ normal(Omega_global, kappa)} (centered;
##' each study's correlation matrix is rejected if not positive definite)
##' rather than a free blend, and the scales use a non-centered log-normal
##' hierarchy, \code{log(tau[i]) ~ normal(ltaum, lsig)}. This
##' parametrization (centered deviations + log-normal tau) was chosen over a
##' non-centered one after an ablation: it removed most divergences and was
##' 1.6-2.7x faster with identical accuracy, except that with weak data
##' (small studies, large \code{kappa}) it can still show some divergences,
##' while mixing better. \code{prior = "rhs"} (the
##' earlier default, kept for comparison/compatibility) is the
##' \emph{regularized horseshoe} of Piironen and Vehtari (2017,
##' \doi{10.1214/17-EJS1337SI}) in a non-centered parametrization, with the
##' regression coefficients also non-centered, blended with each study's own
##' free local correlation via \code{rho}. \code{prior = "horseshoe"} is the
##' earlier ordinary (centered) horseshoe model, kept as a further
##' alternative.
##'
##' \strong{Why \code{"rhs_kappa"} is the default.} The \code{"rhs"} model's
##' \code{rho}-blend has a genuine non-identifiability: because
##' \code{Omega_local[i]} is a free per-study correlation matrix, for any
##' \code{rho} and \code{Omega_global} there is a valid \code{Omega_local[i]}
##' that reproduces study i's data exactly, at zero prior cost. This creates a
##' second posterior mode where \code{Omega_global} collapses toward zero
##' (favoured by the horseshoe) with the true correlation "absorbed"
##' independently into each \code{Omega_local[i]} and \code{rho} left
##' unidentified near its prior mean. This fits studies well but affects
##' exactly what \code{mvn_generate_AP()} uses to simulate an
##' \emph{unobserved} cohort (see \code{mvn_extract_hyperparams()}).
##' \code{"rhs_kappa"} recovered \code{Omega_global} 3-4x more accurately and
##' sampled cleanly every time, while \code{"rhs"} was badly unstable.
##'
##' All three priors need valid initial values (Stan's random defaults fail
##' for \code{"rhs"} and \code{"rhs_kappa"}, whose \code{Omega_global} is
##' built by hand and rejected if not positive definite) and this function
##' supplies them.
##'
##' \strong{Severe-instability check (all priors).} Every returned fit is
##' passed through \code{mvn_diagnose()}; the result is stored in
##' \code{attr(fit, "diagnose")}, and a warning is issued if sampling looks
##' severely unstable (\eqn{\hat R>1.5}, ESS\eqn{{}<5}, or more than 20\% of
##' draws divergent). This is deliberately independent of, and looser than,
##' the \code{"rhs"}-only start check below: that check can only ever compare
##' chains to \emph{each other}, and was found
##' \code{"rhs_kappa"} validation above) to pass silently when two chains
##' agreed with each other while \emph{both} were badly mixing.
##'
##' \strong{Start check (\code{prior = "rhs"} only).} With \code{chains >= 2},
##' \code{check_starts = TRUE} and no user \code{init}, all chains but the last
##' start from data-informed values (pooled within-study residual correlation)
##' while the last chain starts from a deliberately sparse point. If the true
##' correlations are dense, a chain started sparse can fall into a spurious
##' all-zero mode (with reassuring \code{Rhat} if \emph{all} chains do so). If
##' the last chain's posterior-mean \code{Omega_global} differs from the
##' others' by more than \code{start_gap_tol}, and \code{auto_refit = TRUE}
##' (the default), the model is automatically resampled with every chain
##' started from the informed point and a warning names both the original
##' disagreement and the fact that a refit happened (this refit has no alarm
##' chain of its own, so it is not itself protected against the same failure:
##' see \code{auto_refit}). With \code{auto_refit = FALSE} the original
##' (disagreeing) fit is returned with a warning instead. The check's outcome
##' is stored in the \code{"start_check"} attribute of the returned fit either
##' way.
##'
##' @title mvn_infer_mlm_sparse
##' @param Y multivariate response data - each line sign/symptom values for a
##'   patient
##' @param X multivariable explanatory data - each line a set of predictors for
##'   a patient
##' @param study vector of which study each record belongs to
##' @param betaM_prior_sd prior for Betas
##' @param betaS_prior_sd prior for Betas
##' @param tauM_prior_sd prior for tau: SD of the normal prior on the mean
##'   scale (for \code{"rhs_kappa"}: on the mean of log(tau))
##' @param tauS_prior_sd prior for tau: SD of the half-normal prior on the
##'   between-study scale SD (for \code{"rhs_kappa"}: on the SD of log(tau))
##' @param lkj_local_prior_scale prior for local (per-study) correlation
##' @param rhoA beta parameter for rho
##' @param rhoB beta parameter for rho
##' @param prior \code{"rhs_kappa"} (deviation-penalty model, the default:
##'   see Details for why), \code{"rhs"} (regularized horseshoe with the
##'   earlier rho-blend local layer, kept for comparison), or
##'   \code{"horseshoe"} (ordinary horseshoe, the earliest model)
##' @param p0 (\code{"rhs"}/\code{"rhs_kappa"}) prior guess for the number of
##'   non-zero correlations among the \code{choose(ncol(Y),2)} pairs, strictly
##'   between 0 and that number. Sets the scale of the global shrinkage
##'   parameter, \code{p0 / (D - p0) / sqrt(nrow(Y))}, following Piironen and
##'   Vehtari's recipe with \code{1/sqrt(n)} standing in for their
##'   \code{sigma/sqrt(n)} (a heuristic transplant). Default: 10\% of the
##'   pairs, rounded up (5 for 10 variates). Results were insensitive to 5 vs
##'   30 on clustered truth.
##' @param slab_scale (\code{"rhs"}/\code{"rhs_kappa"}) scale of the
##'   Student-t slab that regularizes large correlations, default 0.5
##'   (correlations lie in (-1,1))
##' @param slab_df (\code{"rhs"}/\code{"rhs_kappa"}) slab degrees of freedom,
##'   default 4
##' @param kappa_prior_scale (\code{"rhs_kappa"} only) scale of the
##'   half-normal prior on \code{kappa}, the shared strength of each study's
##'   correlation deviation from \code{Omega_global}. Default 2 (used
##'   throughout validation; not otherwise tuned).
##' @param check_starts (\code{"rhs"} only) use the last chain as an alarm
##'   chain, see Details. Ignored (with no check) if \code{chains < 2} or
##'   \code{init} is supplied. Not used for \code{"rhs_kappa"}, which did not
##'   need it in testing (see Details).
##' @param start_gap_tol correlation-scale threshold for the start check
##'   warning
##' @param auto_refit (\code{"rhs"} only, and only when \code{check_starts}
##'   fires) automatically resample with informed-only starts (roughly
##'   doubling the wall time of a failing fit) rather than just warning and
##'   returning the disagreeing fit. Default \code{TRUE}.
##' @param init optional Stan \code{init} argument. Default \code{NULL}:
##'   informed starts for \code{"rhs"}/\code{"rhs_kappa"}, Stan's random
##'   inits for \code{"horseshoe"}.
##' @param iter iterations for MCMC, default
##' @param cores number of cores to use
##' @param chains number of chains to use
##' @param ...
##' @return a Stan sample object, always carrying an attribute
##'   \code{"diagnose"} (the \code{mvn_diagnose()} result for this fit; see
##'   Details) and, for \code{"rhs"} with a start check, also an attribute
##'   \code{"start_check"}: a list with \code{gap}, \code{tol},
##'   \code{alarm_chain}, \code{agree} and \code{refit} (whether an automatic
##'   refit happened)
##' @author Pete Dodd
##' @references Piironen J, Vehtari A (2017). Sparsity information and
##'   regularization in the horseshoe and other shrinkage priors. Electronic
##'   Journal of Statistics 11(2):5018-5051. \doi{10.1214/17-EJS1337SI}
##' @export
##' @import rstan
mvn_infer_mlm_sparse <- function(Y, X, study,
                                 betaM_prior_sd = 1, # prior for Betas
                                 betaS_prior_sd = 0.5, # prior for Betas
                                 tauM_prior_sd = 1, # prior for tau
                                 tauS_prior_sd = 0.5, # prior for tau
                                 lkj_local_prior_scale = 3, # prior for local cor
                                 rhoA = 2, # beta parameter for rho
                                 rhoB = 2, # beta parameter for rho
                                 prior = c("rhs_kappa", "rhs", "horseshoe"),
                                 p0 = NULL, slab_scale = 0.5, slab_df = 4,
                                 kappa_prior_scale = 2,
                                 check_starts = TRUE, start_gap_tol = 0.1,
                                 auto_refit = TRUE,
                                 init = NULL,
                                 iter = 2e3, cores = 4, chains = 4, ...) {
  prior <- match.arg(prior)
  ## prepare data (records sorted by study, see .sort_by_study)
  srt <- .sort_by_study(Y, X, study)
  Y <- srt$Y
  X <- srt$X
  study <- srt$study
  if (ncol(Y) < 2) stop("the sparse models need at least 2 variates (columns of Y)")
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
    tauM_prior_sd = tauM_prior_sd, # prior for tau
    tauS_prior_sd = tauS_prior_sd, # prior for tau
    lkj_local_prior_scale = lkj_local_prior_scale, # prior for local cor
    rhoA = rhoA, # beta parameter for rho
    rhoB = rhoB # beta parameter for rho
  )

  if (prior == "horseshoe") {
    args <- list(
      object = stanmodels$mvn_inferH_sparse, data = shdata,
      chains = chains, cores = cores, iter = iter
    )
    if (!is.null(init)) args$init <- init
    fit <- do.call(rstan::sampling, c(args, list(...)))
    return(.attach_stability_check(fit, ""))
  }

  if (prior == "rhs_kappa") {
    DR <- choose(ncol(Y), 2)
    if (is.null(p0)) p0 <- .default_p0(ncol(Y))
    if (!(is.numeric(p0) && length(p0) == 1 && p0 > 0 && p0 < DR)) {
      stop("p0 must be a single number with 0 < p0 < choose(ncol(Y), 2) = ", DR)
    }
    if (!(slab_scale > 0 && slab_df > 0)) stop("slab_scale and slab_df must be positive")
    if (!(is.numeric(kappa_prior_scale) && length(kappa_prior_scale) == 1 && kappa_prior_scale > 0)) {
      stop("kappa_prior_scale must be a single positive number")
    }
    ## no rhoA/rhoB/lkj_local_prior_scale for this model: no rho-blend or per-study
    ## LKJ prior to parametrise (see the deviation-penalty design in Details)
    shdata$rhoA <- NULL
    shdata$rhoB <- NULL
    shdata$lkj_local_prior_scale <- NULL
    shdata$p0 <- p0
    shdata$slab_scale <- slab_scale
    shdata$slab_df <- slab_df
    shdata$kappa_prior_scale <- kappa_prior_scale
    if (is.null(init)) init <- .rhs_kappa_inits(Y, X, study, chains, slab_scale)
    fit <- rstan::sampling(stanmodels$mvn_inferH_sparse_kappa,
      data = shdata,
      chains = chains,
      cores = cores,
      iter = iter,
      init = init, ...
    )
    return(.attach_stability_check(fit, ""))
  }

  ## regularized horseshoe
  DR <- choose(ncol(Y), 2)
  if (is.null(p0)) p0 <- .default_p0(ncol(Y))
  if (!(is.numeric(p0) && length(p0) == 1 && p0 > 0 && p0 < DR)) {
    stop("p0 must be a single number with 0 < p0 < choose(ncol(Y), 2) = ", DR)
  }
  if (!(slab_scale > 0 && slab_df > 0)) stop("slab_scale and slab_df must be positive")
  shdata$p0 <- p0
  shdata$slab_scale <- slab_scale
  shdata$slab_df <- slab_df

  do_check <- is.null(init) && check_starts && chains >= 2
  if (is.null(init)) {
    init <- .rhs_inits(Y, X, study, chains, slab_scale, alarm = do_check)
  }
  fit <- rstan::sampling(stanmodels$mvn_inferH_sparse_rhs,
    data = shdata,
    chains = chains,
    cores = cores,
    iter = iter,
    init = init, ...
  )
  if (do_check) {
    gap <- .start_gap(fit, alarm_chain = chains)
    agree <- is.finite(gap) && gap <= start_gap_tol
    if (agree) {
      attr(fit, "start_check") <- list(
        gap = gap, tol = start_gap_tol, alarm_chain = chains, agree = TRUE, refit = FALSE
      )
    } else if (auto_refit) {
      warning(
        "start check failed: chain ", chains, " (started from a sparse point) and the ",
        "data-informed chains disagreed about Omega_global by up to ", signif(gap, 3),
        " (tolerance ", start_gap_tol, "). The posterior may have a spurious all-zero ",
        "mode. Auto-refitting with informed-only starts (auto_refit = TRUE): this refit ",
        "has no alarm chain of its own, so the same failure would not be caught again; ",
        "compare with prior = \"horseshoe\" or mvn_infer_mlm() if in doubt.",
        call. = FALSE
      )
      init2 <- .rhs_inits(Y, X, study, chains, slab_scale, alarm = FALSE)
      fit <- rstan::sampling(stanmodels$mvn_inferH_sparse_rhs,
        data = shdata,
        chains = chains,
        cores = cores,
        iter = iter,
        init = init2, ...
      )
      attr(fit, "start_check") <- list(
        gap = gap, tol = start_gap_tol, alarm_chain = chains, agree = FALSE, refit = TRUE
      )
    } else {
      attr(fit, "start_check") <- list(
        gap = gap, tol = start_gap_tol, alarm_chain = chains, agree = FALSE, refit = FALSE
      )
      warning(
        "start check failed: chain ", chains, " (started from a sparse point) and the ",
        "data-informed chains disagree about Omega_global by up to ", signif(gap, 3),
        " (tolerance ", start_gap_tol, "). The posterior may have a spurious all-zero ",
        "mode; treat this fit (with its disagreeing chains) as unreliable. Refit with ",
        "check_starts = FALSE for informed starts only, or auto_refit = TRUE (the default) ",
        "to do that automatically, and compare with prior = \"horseshoe\" or mvn_infer_mlm().",
        call. = FALSE
      )
    }
  }
  .attach_stability_check(fit, "")
}
