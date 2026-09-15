## tests for the R-side simulators in R/utilities.R (no Stan involved)

test_that("mvn_simulate produces correctly-shaped output and recovers Beta/Sigma", {
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

test_that("mvn_sample_study produces correctly-shaped output", {
  set.seed(102)
  skip_if_not_installed("trialr")
  Ncovars <- 2
  NV <- 3
  Npats <- 50
  X <- cbind(1, runif(Npats))
  OmegaG <- trialr::rlkjcorr(1, NV, 2)

  out <- mvn_sample_study(
    X = X,
    betag = matrix(1, Ncovars, NV), sigb = matrix(0.01, Ncovars, NV),
    rho = 0.5, taug = rep(1, NV), sigt = rep(0.1, NV),
    OmegaG = OmegaG, lkj_local = 3
  )

  expect_equal(dim(out), c(Npats, NV))
})

test_that("mvn_simulate_studies produces the documented data-frame structure", {
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
    rhoA = 2, rhoB = 2, taug = rep(1, NV), sigt = rep(0.2, NV),
    lkj_local = 3, lkj_global = 2
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), Nstudies * Npats)
  expect_true(all(c("studyno", "obsno") %in% names(out)))
  expect_equal(sort(unique(out$studyno)), 1:Nstudies)
  expect_equal(as.vector(table(out$studyno)), rep(Npats, Nstudies))
  ## obsno should restart at 1 within each study
  expect_equal(out$obsno[out$studyno == 1], 1:Npats)
})
