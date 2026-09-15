## tests for mvn_diagnose() (R/diagnostics.R)

test_that("mvn_diagnose returns a well-formed object and correctly handles structural NaN Rhat/ESS", {
  skip_on_cran()
  set.seed(401)
  X <- cbind(1, rnorm(300))
  Y <- mvn_simulate(
    X, matrix(c(0, 0, 0, 1, -1, 2), 2, 3, byrow = TRUE), diag(1, 3)
  )
  # decent iter/chains so this is a genuinely well-behaved fit
  fit <- mvn_infer(
    Y = Y, X = X, iter = 1000, chains = 2, cores = 1, seed = 1, refresh = 0
  )

  ## sanity check on the fixture itself: L_Omega's structural zeros
  ## mvn_diagnose() must not treat that as a failure
  s <- rstan::summary(fit)$summary
  expect_true(any(is.na(s[, "Rhat"])))

  d <- mvn_diagnose(fit)

  expect_s3_class(d, "muvamere_diagnostics")
  expect_named(d, c(
    "rhat_max", "ess_min", "ess_min_ratio", "n_divergent",
    "n_max_treedepth", "bfmi_min", "n_draws", "ok"
  ))
  expect_true(is.finite(d$rhat_max)) # NOT NaN, despite the fixture having NA Rhat entries
  expect_true(is.finite(d$ess_min))
  expect_equal(d$n_draws, (1000 - 500) * 2)
  expect_true(d$ok) # this fit should genuinely be fine

  ## print method runs without error and mentions "OK"
  expect_output(print(d), "OK")
})


test_that("mvn_diagnose's ok flag responds to its thresholds", {
  skip_on_cran()
  set.seed(402)
  X <- cbind(1, rnorm(200))
  Y <- mvn_simulate(X, matrix(0, 2, 2), diag(1, 2))
  fit <- suppressWarnings(
    mvn_infer(
      Y = Y, X = X, iter = 400, chains = 2, cores = 1, seed = 1, refresh = 0
    )
  )

  ## an impossibly strict Rhat threshold must force ok = FALSE, deterministically,
  ## regardless of how well the actual sampler behaved
  diag_strict <- mvn_diagnose(fit, rhat_threshold = -1)
  expect_false(diag_strict$ok)

  ## a trivially loose set of thresholds should not report problems for a
  ## fit with zero divergences/treedepth issues
  diag_loose <- mvn_diagnose(
    fit,
    rhat_threshold = 100, ess_ratio_threshold = 0, bfmi_threshold = 0
  )
  expect_equal(diag_loose$n_divergent, 0)
  expect_true(diag_loose$ok)
})
