## tests for mvn_infer() (single-study Stan model): parameter recovery,
## the Sigma generated quantity, and the Z-prediction pathway. These run
## real (short) MCMC fits, so they're slower than the pure-R tests but still
## fast enough (a few seconds each) for routine use.

test_that("mvn_infer recovers Beta and exposes Sigma as a generated quantity", {
  skip_on_cran()
  set.seed(201)
  Nobs <- 300
  NP <- 2
  NV <- 3
  X <- cbind(1, rnorm(Nobs))
  Beta_true <- matrix(c(0, 0, 0, 1, -1, 2), nrow = NP, ncol = NV, byrow = TRUE)
  Sigma_true <- diag(1, NV)
  Sigma_true[1, 2] <- Sigma_true[2, 1] <- 0.5
  Y <- mvn_simulate(X, Beta_true, Sigma_true)

  fit <- mvn_infer(
    Y = Y, X = X, iter = 800, chains = 2, cores = 1, seed = 1, refresh = 0
  )

  Beta_hat <- matrix(
    rstan::summary(fit, pars = "Beta")$summary[, "mean"],
    nrow = NP, ncol = NV, byrow = TRUE
  )
  expect_equal(Beta_hat, Beta_true, tolerance = 0.35)

  Sigma_hat <- matrix(
    rstan::summary(fit, pars = "Sigma")$summary[, "mean"], NV, NV
  )
  expect_equal(Sigma_hat, Sigma_true, tolerance = 0.35)
})


test_that("mvn_infer rejects a Z with the wrong number of columns", {
  set.seed(202)
  Nobs <- 50
  NP <- 2
  NV <- 3
  X <- cbind(1, rnorm(Nobs))
  Y <- mvn_simulate(X, matrix(0, NP, NV), diag(1, NV))

  expect_error(
    mvn_infer(
      Y = Y, X = X, Z = matrix(1, 5, NV),
      iter = 10, chains = 1, refresh = 0
    ),
    "same number of columns"
  )
})


test_that("mvn_infer's Z-prediction pathway works when NP != NV, and mvn_extract_predictions returns the right shape", {
  skip_on_cran()
  set.seed(203)
  Nobs <- 300
  NP <- 2
  NV <- 5 # deliberately != NP -- this exact case used to be silently broken
  X <- cbind(1, rnorm(Nobs))
  Beta_true <- matrix(
    c(0, 0, 0, 0, 0, 2, -1, 1, 0.5, -0.5),
    nrow = NP, ncol = NV, byrow = TRUE
  )
  Sigma_true <- diag(1, NV)
  Y <- mvn_simulate(X, Beta_true, Sigma_true)

  Znew <- cbind(1, c(-1, 0, 1))
  fit <- mvn_infer(
    Y = Y, X = X, Z = Znew,
    iter = 800, chains = 2, cores = 1, seed = 1, refresh = 0
  )

  preds <- mvn_extract_predictions(fit)
  expect_equal(dim(preds), c(800, nrow(Znew), NV)) # [ndraws, NewObs, NV]

  pred_mean <- apply(preds, c(2, 3), mean)
  true_mean <- Znew %*% Beta_true
  expect_equal(pred_mean, true_mean, tolerance = 0.5, ignore_attr = TRUE)
})

test_that("mvn_extract_predictions errors informatively when Z was not supplied", {
  skip_on_cran()
  set.seed(204)
  Nobs <- 50
  NP <- 2
  NV <- 2
  X <- cbind(1, rnorm(Nobs))
  Y <- mvn_simulate(X, matrix(0, NP, NV), diag(1, NV))
  ## small iter here is expected to produce sampler warnings (Rhat/ESS) that
  ## are irrelevant to what this test checks (the no-Z error message)
  fit <- suppressWarnings(mvn_infer(
    Y = Y, X = X,
    iter = 200, chains = 1, refresh = 0
  ))

  expect_error(mvn_extract_predictions(fit), "non-NULL Z")
})
