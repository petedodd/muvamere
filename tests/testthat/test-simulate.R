## tests for the R-side simulators in R/utilities.R (no Stan involved)

test_that("mvn_simulate output is well-shaped and recovers Beta/Sigma", {
  set.seed(101)
  Nobs <- 4000
  NP <- 2
  NV <- 3
  X <- cbind(1, rnorm(Nobs))
  Beta <- matrix(c(0, 0, 0, 2, -1, 0.5), nrow = NP, ncol = NV, byrow = TRUE)
  Sigma <- diag(1, NV)
  Sigma[1, 2] <- Sigma[2, 1] <- 0.4

  Y <- mvn_simulate(X, Beta, Sigma)

  expect_equal(dim(Y), c(Nobs, NV))
  ## with Nobs=4000 the empirical mean/cov should be close to the true values
  ## generous tolerance since this is a stochastic check, not exact
  fit_lm_means <- colMeans(Y - X %*% Beta) # residual means should be ~0
  expect_true(all(abs(fit_lm_means) < 0.1))
  expect_equal(cov(Y - X %*% Beta), Sigma, tolerance = 0.15)
})

test_that("mvn_sample_study_kappa gives log-normal scales, right shape", {
  set.seed(102)
  NV <- 3
  Npats <- 4000
  X <- cbind(1, runif(Npats))
  ## kappa = 0 and sigb = 0: the only randomness in scale is log-normal tau,
  ## and with lsig = 0 every tau is exactly exp(ltaum)
  out <- mvn_sample_study_kappa(X,
    betag = matrix(1, 2, NV), sigb = matrix(0, 2, NV), kappa = 0,
    ltaum = log(c(0.5, 1, 2)), lsig = rep(0, NV), OmegaG = diag(NV)
  )
  expect_equal(dim(out), c(Npats, NV))
  expect_equal(apply(out - X %*% matrix(1, 2, NV), 2, sd), c(0.5, 1, 2),
    tolerance = 0.05
  )
})

test_that("mvn_simulate_studies returns the documented data-frame structure", {
  set.seed(103)
  skip_if_not_installed("trialr")
  Ncovars <- 2
  NV <- 2
  Npats <- 20
  Nstudies <- 4
  Xlist <- replicate(Nstudies, cbind(1, runif(Npats)), simplify = FALSE)

  out <- mvn_simulate_studies(
    Xlist,
    betag = matrix(1, Ncovars, NV), sigb = matrix(0.05, Ncovars, NV),
    kappa = 0.05, ltaum = rep(0, NV), lsig = rep(0.1, NV), lkj_global = 2
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), Nstudies * Npats)
  expect_true(all(c("studyno", "obsno") %in% names(out)))
  expect_equal(sort(unique(out$studyno)), 1:Nstudies)
  expect_equal(as.vector(table(out$studyno)), rep(Npats, Nstudies))
  ## obsno should restart at 1 within each study
  expect_equal(out$obsno[out$studyno == 1], 1:Npats)
  ## the drawn global correlation is returned as an attribute
  expect_equal(dim(attr(out, "OmegaG")), c(NV, NV))
})
