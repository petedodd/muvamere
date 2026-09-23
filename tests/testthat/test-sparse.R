## tests for mvn_infer_mlm_sparse(): the deviation-penalty (rhs_kappa,
## default), regularized-horseshoe (rhs) and ordinary-horseshoe models, the
## informed initial values, the start-agreement check (rhs only) and the
## general severe-instability check (all priors).

## simulate multi-study data with a KNOWN global correlation matrix OG
## (mvn_simulate_studies() draws its own random global correlation, so it
## cannot be used to test recovery of a specific one)
sim_known_global <- function(OG, S = 3, Np = 150, seed = 1) {
  set.seed(seed)
  NV <- ncol(OG)
  Xlist <- replicate(S, cbind(1, runif(Np)), simplify = FALSE)
  Ys <- lapply(Xlist, function(X) {
    mvn_sample_study(
      X,
      matrix(1, 2, NV),
      matrix(0.05, 2, NV),
      0.9,
      rep(1, NV),
      rep(0.1, NV),
      OG,
      3
    )
  })
  list(
    Y = do.call(rbind, Ys),
    X = do.call(rbind, Xlist),
    study = rep(seq_len(S), each = Np)
  )
}


test_that(".default_p0 validly guesses the number of non-zero correlations", {
  expect_equal(
    muvamere:::.default_p0(10), 5
  ) # 10% of 45 pairs, rounded up
  for (nv in 2:12) {
    p0 <- muvamere:::.default_p0(nv)
    expect_true(
      p0 > 0 && p0 < choose(nv, 2)
    ) # Stan requires 0 < p0 < choose(NV,2)
  }
})

test_that(".upper_idx matches Stan's row-major upper-triangle numbering", {
  idx <- muvamere:::.upper_idx(4)
  expect_equal(nrow(idx), 6)
  expect_equal(unname(idx[1:4, ]), rbind(c(1, 2), c(1, 3), c(1, 4), c(2, 3)))
  expect_true(all(idx[, 1] < idx[, 2]))
})

test_that(".rhs_inits gives PD informed starts and a sparse alarm chain", {
  set.seed(1)
  OG <- matrix(0.4, 4, 4); diag(OG) <- 1
  d <- sim_known_global(OG, S = 2, Np = 100)
  ini <- muvamere:::.rhs_inits(
    d$Y, d$X, d$study,
    chains = 3, slab_scale = 0.5, alarm = TRUE
  )
  expect_length(ini, 3)
  ## informed chains: rebuild Omega_global from (z, T, lam = 1, caux = 1)
  ## exactly as the Stan model does
  rebuild <- function(i) {
    T0 <- i$T; c0 <- 0.5
    lt <- c0 / sqrt(c0^2 + T0^2)
    idx <- muvamere:::.upper_idx(4)
    O <- diag(4)
    for (k in seq_len(nrow(idx))) {
      O[idx[k, 1], idx[k, 2]] <- O[idx[k, 2], idx[k, 1]] <- i$zg[k] * T0 * lt
    }
    O
  }
  for (ch in 1:2) {
    expect_equal(ini[[ch]]$T, 0.3)
    O <- rebuild(ini[[ch]])
    expect_true(
      min(eigen(O)$values) > 0
    ) # positive definite, or Stan would reject the start
    expect_true(
      mean(O[upper.tri(O)]) > 0.2
    ) # carries the (dense, positive) sample correlation
  }
  expect_equal(ini[[3]]$T, 0.01) # alarm chain = sparse start
  ## without an alarm every chain is informed
  ini2 <- muvamere:::.rhs_inits(
    d$Y, d$X, d$study,
    chains = 2, slab_scale = 0.5, alarm = FALSE
  )
  expect_equal(
    vapply(ini2, function(i) i$T, numeric(1)), c(0.3, 0.3)
  )
})

test_that("mvn_infer_mlm_sparse validates p0 and #variates before sampling", {
  d <- sim_known_global(diag(3), S = 2, Np = 30)
  expect_error(
    mvn_infer_mlm_sparse(d$Y, d$X, d$study, p0 = 3), "p0 must be"
  )
  expect_error(
    mvn_infer_mlm_sparse(d$Y, d$X, d$study, p0 = 0), "p0 must be"
  )
  expect_error(
    mvn_infer_mlm_sparse(d$Y[, 1, drop = FALSE], d$X, d$study),
    "at least 2 variates"
  )
  expect_error(
    mvn_infer_mlm_sparse(d$Y, d$X, d$study, prior = "nope"), "should be one of"
  )
})

test_that("prior = 'rhs' recovers a DENSE global correlation (not all-zero)", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  OG <- matrix(0.4, 4, 4); diag(OG) <- 1
  d <- sim_known_global(OG, S = 3, Np = 150, seed = 11)
  fit <- suppressWarnings(mvn_infer_mlm_sparse(
    d$Y, d$X, d$study,
    prior = "rhs", iter = 600, chains = 2, cores = 1, refresh = 0, seed = 1
  ))
  expect_equal(fit@par_dims$Omega_global, c(4, 4))
  expect_equal(fit@par_dims$Omega_local, c(3, 4, 4))
  expect_equal(fit@par_dims$Betas, c(3, 2, 4))
  hyper <- mvn_extract_hyperparams(fit)
  Og <- hyper$OmegaG
  expect_equal(diag(Og), rep(1, 4), tolerance = 1e-6)
  expect_true(min(eigen(Og)$values) > 0)
  ## truth is 0.4 everywhere; the failure this guards against estimates ~0
  expect_true(mean(Og[upper.tri(Og)]) > 0.25)
  expect_true(
    hyper$rho > 0.5
  ) # rho is data-informed (truth 0.9), not stuck at its prior mean ~0.5
  chk <- attr(fit, "start_check")
  expect_true(chk$agree)
  expect_equal(chk$alarm_chain, 2)
})

test_that("prior = 'rhs' shrinks a sparse global correlation toward zero", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  OG <- diag(4)
  OG[1, 2] <- OG[2, 1] <- 0.7
  d <- sim_known_global(OG, S = 3, Np = 150, seed = 12)
  fit <- suppressWarnings(mvn_infer_mlm_sparse(
    d$Y, d$X, d$study,
    prior = "rhs", iter = 600, chains = 2, cores = 1, refresh = 0, seed = 1
  ))
  Og <- mvn_extract_hyperparams(fit)$OmegaG
  expect_true(abs(Og[1, 2] - 0.7) < 0.15)
  zeros <- Og[
    upper.tri(Og)
  ][-1] # every pair except (1,2), whose column-major position is first
  expect_true(all(
    abs(zeros) < 0.15
  ))
})

test_that("the start check raises a warning when the alarm chain disagrees", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  OG <- diag(3)
  d <- sim_known_global(OG, S = 2, Np = 60, seed = 13)
  ## tolerance 0 => any difference at all counts as disagreement: exercises the
  ## warning path
  warns <- character()
  fit <- withCallingHandlers(
    mvn_infer_mlm_sparse(d$Y, d$X, d$study, prior = "rhs",
      iter = 200, chains = 2, cores = 1,
      refresh = 0, seed = 1, start_gap_tol = 0
    ),
    ## collect (and silence) all warnings, incl. rstan's short-run ESS ones
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("start check failed", warns)))
  ## auto_refit = TRUE is the default
  expect_true(any(grepl("Auto-refitting", warns)))
  chk <- attr(fit, "start_check")
  expect_false(chk$agree)
  ## the returned fit is the informed-only refit, not the disagreeing one
  expect_true(chk$refit)
  ## and no check with a single chain, a user init, or check_starts = FALSE
  f1 <- suppressWarnings(mvn_infer_mlm_sparse(d$Y, d$X, d$study, prior = "rhs",
    iter = 200, chains = 1, cores = 1, refresh = 0, seed = 1
  ))
  expect_null(attr(f1, "start_check"))
  f2 <- suppressWarnings(mvn_infer_mlm_sparse(d$Y, d$X, d$study, prior = "rhs",
    iter = 200, chains = 2, cores = 1,
    refresh = 0, seed = 1, check_starts = FALSE
  ))
  expect_null(attr(f2, "start_check"))
})

test_that("auto_refit = FALSE returns the disagreeing fit, no resampling", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  OG <- diag(3)
  d <- sim_known_global(OG, S = 2, Np = 60, seed = 13)
  warns <- character()
  fit <- withCallingHandlers(
    mvn_infer_mlm_sparse(d$Y, d$X, d$study, prior = "rhs",
      iter = 200, chains = 2, cores = 1,
      refresh = 0, seed = 1, start_gap_tol = 0, auto_refit = FALSE
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("start check failed", warns)))
  expect_false(any(grepl("Auto-refitting", warns)))
  chk <- attr(fit, "start_check")
  expect_false(chk$agree)
  expect_false(chk$refit)
})

test_that("prior = 'horseshoe' still fits the ordinary-horseshoe model", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  d <- sim_known_global(diag(3), S = 2, Np = 60, seed = 14)
  fit <- suppressWarnings(mvn_infer_mlm_sparse(d$Y, d$X, d$study,
    prior = "horseshoe", iter = 200, chains = 1, cores = 1,
    refresh = 0, seed = 1
  ))
  expect_equal(
    fit@par_dims$Omega_global, c(3, 3)
  ) # a corr_matrix parameter here
  expect_false("zg" %in% names(fit@par_dims)) # not the non-centered model
})

test_that("prior = 'rhs_kappa' fits and recovers a known correlation", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  OG <- diag(4); OG[1, 2] <- OG[2, 1] <- 0.6
  d <- sim_known_global(OG, S = 3, Np = 150, seed = 21)
  fit <- suppressWarnings(mvn_infer_mlm_sparse(d$Y, d$X, d$study,
    prior = "rhs_kappa", iter = 400, chains = 2, cores = 1, refresh = 0,
    seed = 1
  ))
  expect_equal(fit@par_dims$Omega_global, c(4, 4))
  expect_true("kappa" %in% names(fit@par_dims))
  expect_false("rho" %in% names(fit@par_dims)) # no rho-blend in this model
  ## centered per-study correlations and the log-normal tau hierarchy
  expect_equal(fit@par_dims$rr, c(3, 6))
  expect_false("eta" %in% names(fit@par_dims))
  expect_true(all(c("ltaum", "lsig") %in% names(fit@par_dims)))
  ## no free per-study correlation either
  expect_false("Omega_local" %in% names(fit@par_dims))
  hyper <- mvn_extract_hyperparams(fit)
  expect_equal(hyper$model, "kappa")
  expect_true(is.numeric(hyper$kappa) && hyper$kappa >= 0)
  expect_length(hyper$ltaum, 4)
  expect_true(all(hyper$lsig >= 0))
  expect_null(hyper$taug) # truncated-normal tau is the rho models' only
  Og <- hyper$OmegaG
  expect_equal(diag(Og), rep(1, 4), tolerance = 1e-6)
  expect_true(min(eigen(Og, symmetric = TRUE, only.values = TRUE)$values) > 0)
  expect_true(abs(Og[1, 2] - 0.6) < 0.2)
  ## the other 5 true-zero pairs stay small
  expect_true(all(abs(Og[upper.tri(Og)][-1]) < 0.2))
  ## no start_check attribute for this prior (validated separately not to need
  ## one)
  expect_null(attr(fit, "start_check"))
  ## but the general stability check always attaches
  expect_s3_class(attr(fit, "diagnose"), "muvamere_diagnostics")
})

test_that("'rhs_kappa' validates kappa_prior_scale; drops rho/lkj_local_*", {
  d <- sim_known_global(diag(3), S = 2, Np = 30)
  expect_error(
    mvn_infer_mlm_sparse(d$Y, d$X, d$study,
      prior = "rhs_kappa", kappa_prior_scale = 0
    ),
    "kappa_prior_scale"
  )
  expect_error(
    mvn_infer_mlm_sparse(d$Y, d$X, d$study,
      prior = "rhs_kappa", kappa_prior_scale = c(1, 2)
    ),
    "kappa_prior_scale"
  )
})

test_that("mvn_generate_AP dispatches to the kappa generator for rhs_kappa", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  OG <- diag(4); OG[1, 2] <- OG[2, 1] <- 0.6
  d <- sim_known_global(OG, S = 3, Np = 150, seed = 22)
  fit <- suppressWarnings(mvn_infer_mlm_sparse(d$Y, d$X, d$study,
    prior = "rhs_kappa", iter = 400, chains = 2, cores = 1, refresh = 0,
    seed = 1
  ))
  Xtarget <- replicate(2, cbind(1, runif(10)), simplify = FALSE)
  AP <- mvn_generate_AP(fit, Xtarget)
  expect_s3_class(AP, "data.frame")
  expect_equal(nrow(AP), 20)
  expect_true(all(c("studyno", "obsno") %in% names(AP)))
})

test_that("mvn_sample_study_kappa gives valid, well-shaped correlation", {
  set.seed(31)
  OG <- diag(4); OG[1, 2] <- OG[2, 1] <- 0.5
  X <- cbind(1, runif(200))
  Y <- mvn_sample_study_kappa(X, matrix(1, 2, 4), matrix(0.1, 2, 4),
    kappa = 0.05,
    ltaum = rep(0, 4), lsig = rep(0.1, 4), OmegaG = OG
  )
  expect_equal(dim(Y), c(200, 4))
  ## a tiny kappa means the generated study's own empirical correlation should
  ## still look like OG
  r <- cor(Y)
  expect_true(abs(r[1, 2] - 0.5) < 0.3)
})

test_that(".attach_stability_check flags instability despite agreement", {
  ## a fabricated mvn_diagnose()-shaped list standing in for a
  ## badly-mixing-but-agreeing fit
  bad <- structure(
    list(rhat_max = 5, ess_min = 2, ess_min_ratio = 0.002, n_divergent = 400,
      n_draws = 1000, n_max_treedepth = 0, bfmi_min = 0.5, ok = FALSE),
    class = "muvamere_diagnostics"
  )
  good <- structure(
    ## ok = FALSE (mvn_diagnose's strict bar) but NOT severely unstable by the
    ## looser bar used here
    list(rhat_max = 1.06, ess_min = 40, ess_min_ratio = 0.1, n_divergent = 8,
      n_draws = 1000, n_max_treedepth = 0, bfmi_min = 0.5, ok = FALSE),
    class = "muvamere_diagnostics"
  )
  expect_true(muvamere:::.severely_unstable(bad))
  expect_false(muvamere:::.severely_unstable(good))
})

test_that(".rhs_kappa_inits starts study correlations at Omega_global", {
  set.seed(3)
  d <- sim_known_global(diag(4), S = 3, Np = 50)
  ini <- muvamere:::.rhs_kappa_inits(d$Y, d$X, d$study, chains = 2,
    slab_scale = 0.5
  )
  expect_length(ini, 2)
  for (l in ini) {
    expect_null(l$eta)
    expect_equal(dim(l$rr), c(3, 6)) # [studies, pairs], a plain matrix
    ## identical rows = Omega_global's upper triangle, so every study's
    ## correlation matrix starts positive definite
    og <- l$zg * 0.3 * 0.5 / sqrt(0.5^2 + 0.3^2)
    for (i in 1:3) expect_equal(l$rr[i, ], og)
  }
})
