## tests for mvn_infer_mlm() (hierarchical Stan model) and the R/bridge.R
## functions built on top of it.

test_that("mvn_infer_mlm derives Nstudies from the number of distinct studies, not Nrecords", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  set.seed(301)
  Nstudies <- 3
  Npats <- 25 # Nrecords = Nstudies * Npats = 75 -- must NOT show up as par_dims$tau[1]
  NV <- 2
  Ncov <- 2
  Xlist <- replicate(Nstudies, cbind(1, runif(Npats)), simplify = FALSE)

  allsamp <- mvn_simulate_studies(
    Xlist,
    betag = matrix(1, Ncov, NV), sigb = matrix(0.05, Ncov, NV),
    rhoA = 2, rhoB = 2, taug = rep(1, NV), sigt = rep(0.2, NV),
    lkj_local = 3, lkj_global = 2
  )
  expect_equal(
    nrow(allsamp), Nstudies * Npats
  ) # sanity check on the fixture itself

  fit <- suppressWarnings(mvn_infer_mlm(
    Y = as.matrix(allsamp[, 1:NV]), X = do.call(rbind, Xlist),
    study = allsamp$studyno, iter = 300, chains = 1, cores = 1, refresh = 0
  ))

  expect_equal(fit@par_dims$tau, c(Nstudies, NV))
  expect_equal(fit@par_dims$Omega_local, c(Nstudies, NV, NV))
  expect_equal(fit@par_dims$Betas, c(Nstudies, Ncov, NV))
})

test_that("mvn_extract_hyperparams returns correctly-shaped, correctly-reshaped hyperparameters", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  set.seed(302)
  Nstudies <- 3
  Npats <- 25
  NV <- 2
  Ncov <- 2
  Xlist <- replicate(Nstudies, cbind(1, runif(Npats)), simplify = FALSE)
  allsamp <- mvn_simulate_studies(
    Xlist,
    betag = matrix(1, Ncov, NV), sigb = matrix(0.05, Ncov, NV),
    rhoA = 2, rhoB = 2, taug = rep(1, NV), sigt = rep(0.2, NV),
    lkj_local = 3, lkj_global = 2
  )
  fit <- suppressWarnings(mvn_infer_mlm(
    Y = as.matrix(allsamp[, 1:NV]), X = do.call(rbind, Xlist),
    study = allsamp$studyno, iter = 300, chains = 1, cores = 1, refresh = 0
  ))

  hyper <- mvn_extract_hyperparams(fit)

  expect_named(hyper, c("betag", "sigb", "taug", "sigt", "rho", "OmegaG"))
  expect_equal(dim(hyper$betag), c(Ncov, NV))
  expect_equal(dim(hyper$sigb), c(Ncov, NV))
  expect_length(hyper$taug, NV)
  expect_length(hyper$sigt, NV)
  expect_length(hyper$rho, 1)
  expect_true(hyper$rho >= 0 && hyper$rho <= 1)
  expect_equal(dim(hyper$OmegaG), c(NV, NV))
  expect_equal(diag(hyper$OmegaG), rep(1, NV), tolerance = 1e-6) # a correlation matrix
})

test_that("mvn_generate_AP produces output in the same format as mvn_simulate_studies", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  set.seed(303)
  Nstudies <- 3
  Npats <- 25
  NV <- 2
  Ncov <- 2
  Xlist <- replicate(Nstudies, cbind(1, runif(Npats)), simplify = FALSE)
  allsamp <- mvn_simulate_studies(
    Xlist,
    betag = matrix(1, Ncov, NV), sigb = matrix(0.05, Ncov, NV),
    rhoA = 2, rhoB = 2, taug = rep(1, NV), sigt = rep(0.2, NV),
    lkj_local = 3, lkj_global = 2
  )
  fit <- suppressWarnings(mvn_infer_mlm(
    Y = as.matrix(allsamp[, 1:NV]), X = do.call(rbind, Xlist),
    study = allsamp$studyno, iter = 300, chains = 1, cores = 1, refresh = 0
  ))

  ## generate for a different (here: same-shaped) set of target studies
  Xtarget <- replicate(2, cbind(1, runif(10)), simplify = FALSE)
  AP <- mvn_generate_AP(fit, Xtarget, lkj_local = 3)

  expect_s3_class(AP, "data.frame")
  expect_equal(nrow(AP), 20)
  expect_true(all(c("studyno", "obsno") %in% names(AP)))
  expect_equal(sort(unique(AP$studyno)), 1:2)
  expect_equal(as.vector(table(AP$studyno)), c(10, 10))
})


## the hierarchical Stan models evaluate the likelihood study-by-study and need
## records sorted by study; the R wrappers sort via .sort_by_study().

test_that(".sort_by_study sorts records into contiguous study blocks and relabels ids 1..S", {
  Y <- cbind(a = 1:6, b = 11:16)
  X <- cbind(1, c(6, 5, 4, 3, 2, 1))
  study <- c(20, 10, 20, 30, 10, 30) # unsorted, ids not 1..S
  s <- muvamere:::.sort_by_study(Y, X, study)
  expect_equal(s$study, c(1L, 1L, 2L, 2L, 3L, 3L))
  expect_equal(s$Y[, "a"], c(2L, 5L, 1L, 3L, 4L, 6L)) # rows moved with their study
  expect_equal(s$X[, 2], X[order(study), 2]) # X moved consistently with Y
  expect_error(muvamere:::.sort_by_study(Y, X, study[-1]), "same number of records")
})

test_that("mvn_infer_mlm gives the same fit whatever the record order", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  set.seed(302)
  NV <- 2; Ncov <- 2
  Xlist <- list(cbind(1, runif(20)), cbind(1, runif(35)))
  allsamp <- mvn_simulate_studies(
    Xlist,
    betag = matrix(1, Ncov, NV), sigb = matrix(0.05, Ncov, NV),
    rhoA = 2, rhoB = 2, taug = rep(1, NV), sigt = rep(0.2, NV),
    lkj_local = 3, lkj_global = 2
  )
  Y <- as.matrix(allsamp[, 1:NV]); X <- do.call(rbind, Xlist); st <- allsamp$studyno
  ## interleave the studies but keep each study's own record order, so that after
  ## sorting the data are bitwise identical
  perm <- order(ave(seq_along(st), st, FUN = seq_along) + runif(length(st), 0, 0.5))
  expect_false(is.unsorted(perm) == FALSE) # fixture really is reordered
  fit1 <- suppressWarnings(mvn_infer_mlm(Y, X, st, iter = 200, chains = 1, cores = 1, refresh = 0, seed = 5))
  fit2 <- suppressWarnings(mvn_infer_mlm(Y[perm, ], X[perm, ], st[perm], iter = 200, chains = 1, cores = 1, refresh = 0, seed = 5))
  ## identical data after sorting, same seed => identical draws (permuted=FALSE:
  ## the default permuted=TRUE shuffles draws with R's RNG, which would differ)
  expect_equal(
    rstan::extract(fit1, pars = "rho", permuted = FALSE),
    rstan::extract(fit2, pars = "rho", permuted = FALSE)
  )
})
