## tests for mvn_infer_mlm_sparse(), the deviation-penalty (rhs_kappa) model:
## argument checks, informed initial values, recovery, generation, and the
## severe-instability check. (Tests of the removed rho-blend / horseshoe
## models and their start-agreement check were deleted 2026-09-24.)

## simulate multi-study data with a KNOWN global correlation matrix OG
## (mvn_simulate_studies() draws its own random global correlation, so it
## cannot be used to test recovery of a specific one)
sim_known_global <- function(OG, S = 3, Np = 150, seed = 1) {
  set.seed(seed)
  NV <- ncol(OG)
  Xlist <- replicate(S, cbind(1, runif(Np)), simplify = FALSE)
  Ys <- lapply(Xlist, function(X) {
    mvn_sample_study_kappa(X,
      betag = matrix(1, 2, NV), sigb = matrix(0.05, 2, NV), kappa = 0.05,
      ltaum = rep(0, NV), lsig = rep(0.1, NV), OmegaG = OG
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
    mvn_infer_mlm_sparse(d$Y, d$X, d$study, prior = "rhs"), "were removed"
  )
  expect_error(
    mvn_infer_mlm_sparse(d$Y, d$X, d$study, prior = "horseshoe"), "removed"
  )
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
  Og <- hyper$OmegaG
  expect_equal(unname(diag(Og)), rep(1, 4), tolerance = 1e-6)
  expect_true(min(eigen(Og, symmetric = TRUE, only.values = TRUE)$values) > 0)
  expect_true(abs(Og[1, 2] - 0.6) < 0.2)
  ## the other 5 true-zero pairs stay small
  expect_true(all(abs(Og[upper.tri(Og)][-1]) < 0.2))
  ## per-study Omega_s/Sigs are not saved by default
  expect_false(any(grepl("^(Omega_s|Sigs)\\[", names(fit))))
  expect_true(any(grepl("^Omega_global\\[", names(fit))))
  ## the general stability check always attaches
  expect_s3_class(attr(fit, "diagnose"), "muvamere_diagnostics")
})

test_that("mvn_infer_mlm_sparse validates kappa_prior_scale", {
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

test_that("mvn_generate_AP generates new studies from an rhs_kappa fit", {
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

test_that(".severely_unstable flags instability even when chains agree", {
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
    og <- as.vector(l$zg) * 0.3 * 0.5 / sqrt(0.5^2 + 0.3^2)
    for (i in 1:3) expect_equal(l$rr[i, ], og)
  }
})

test_that("mvn_infer_mlm_sparse runs with only 2 variates (1 correlation)", {
  skip_on_cran()
  ## regression: length-1 zg/lam inits were read by rstan as scalars
  d <- sim_known_global(diag(2), S = 2, Np = 40)
  fit <- suppressWarnings(mvn_infer_mlm_sparse(d$Y, d$X, d$study,
    iter = 100, chains = 1, cores = 1, refresh = 0, seed = 1
  ))
  expect_equal(fit@par_dims$zg, 1)
  expect_s3_class(attr(fit, "diagnose"), "muvamere_diagnostics")
})
