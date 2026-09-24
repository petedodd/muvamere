## tests for the mixed continuous + binary model (mvn_infer_mlm_mixed()), the
## binary-emitting simulators, and the data-free generator object.

## mixed data with a known global correlation: variates (c1, b1, c2, b2)
## in a deliberately interleaved user order, binaries b1/b2
sim_mixed <- function(S = 3, Np = 150, seed = 1, r = 0.6) {
  set.seed(seed)
  OG <- diag(4)
  OG[1, 2] <- OG[2, 1] <- r # c1 with b1 (biserial)
  binary <- c(FALSE, TRUE, FALSE, TRUE)
  betag <- matrix(c(1, 0.2, 1, -0.3, 0.5, 0.3, 1, 0.1), 2, 4)
  colnames(betag) <- c("c1", "b1", "c2", "b2")
  Xlist <- replicate(S, cbind(1, runif(Np)), simplify = FALSE)
  Ys <- lapply(Xlist, function(X) {
    mvn_sample_study_kappa(X, betag, matrix(0.05, 2, 4), kappa = 0.05,
      ltaum = c(0, 0), lsig = c(0.1, 0.1), OmegaG = OG, binary = binary)
  })
  list(Y = do.call(rbind, Ys), X = do.call(rbind, Xlist),
    study = rep(seq_len(S), each = Np), OG = OG, binary = binary)
}

test_that("mvn_sample_study_kappa emits 0/1 binaries with latent scale 1", {
  set.seed(5)
  X <- cbind(1, runif(5000))
  betag <- matrix(0, 2, 3, dimnames = list(NULL, c("a", "b", "c")))
  Y <- mvn_sample_study_kappa(X, betag, matrix(0, 2, 3), kappa = 0,
    ltaum = log(2), lsig = 0, OmegaG = diag(3),
    binary = c(FALSE, TRUE, TRUE))
  expect_equal(colnames(Y), c("a", "b", "c"))
  expect_true(all(Y[, 2:3] %in% c(0, 1)))
  ## latent mean 0 => P(1) = 0.5; continuous scale exp(log 2) = 2
  expect_equal(mean(Y[, 2]), 0.5, tolerance = 0.03)
  expect_equal(sd(Y[, 1]), 2, tolerance = 0.05)
  expect_error(mvn_sample_study_kappa(X, betag, matrix(0, 2, 3), 0, 0, 0,
    diag(3), binary = TRUE), "one entry per variate")
})

test_that("mvn_simulate_studies passes binary through", {
  set.seed(6)
  skip_if_not_installed("trialr")
  out <- mvn_simulate_studies(list(cbind(1, runif(20)), cbind(1, runif(30))),
    betag = matrix(0, 2, 3), sigb = matrix(0.1, 2, 3), kappa = 0.05,
    ltaum = 0, lsig = 0.1, lkj_global = 2, binary = c(TRUE, FALSE, TRUE))
  expect_true(all(unlist(out[, c(1, 3)]) %in% c(0, 1)))
  expect_false(all(out[, 2] %in% c(0, 1)))
})

test_that(".binary_cols accepts logical, indices and names", {
  Y <- matrix(0, 2, 3, dimnames = list(NULL, c("a", "b", "c")))
  f <- muvamere:::.binary_cols
  expect_equal(f(c(FALSE, TRUE, TRUE), Y), c(FALSE, TRUE, TRUE))
  expect_equal(f(2:3, Y), c(FALSE, TRUE, TRUE))
  expect_equal(f(c("c", "b"), Y), c(FALSE, TRUE, TRUE))
  expect_error(f("z", Y), "not in Y")
  expect_error(f(TRUE, Y), "ncol")
})

test_that("mvn_infer_mlm_mixed rejects bad input before sampling", {
  d <- sim_mixed(S = 2, Np = 40)
  expect_error(mvn_infer_mlm_mixed(d$Y, d$X, d$study,
    binary = rep(FALSE, 4)), "no binary columns")
  Y <- d$Y; Y[1, 2] <- 2
  expect_error(mvn_infer_mlm_mixed(Y, d$X, d$study, binary = d$binary),
    "only 0, 1")
  X <- d$X; X[1, 2] <- NA
  expect_error(mvn_infer_mlm_mixed(d$Y, X, d$study, binary = d$binary),
    "X must not contain missing")
  ## guard rail: b2 entirely missing in study 1 -> refused
  Y <- d$Y; Y[d$study == 1, 4] <- NA
  expect_error(mvn_infer_mlm_mixed(Y, d$X, d$study, binary = d$binary),
    "too few jointly observed")
  ## ... and so is a variable pair observed together too rarely
  ## (study 2: c1 only on even rows, b1 only on odd rows -> never jointly)
  i2 <- which(d$study == 2)
  Y <- d$Y
  Y[i2[i2 %% 2 == 0], "c1"] <- NA
  Y[i2[i2 %% 2 == 1], "b1"] <- NA
  expect_error(mvn_infer_mlm_mixed(Y, d$X, d$study, binary = d$binary),
    "c1 & b1")
  ## the continuous wrapper refuses missing values
  Y <- d$Y; Y[1, 1] <- NA
  expect_error(mvn_infer_mlm_sparse(Y[, c(1, 3)], d$X, d$study),
    "mvn_infer_mlm_mixed")
})

## one shared small fit for the tests below (fitting is the slow part)
fit_mixed <- NULL
get_fit <- function() {
  if (is.null(fit_mixed)) {
    d <- sim_mixed(S = 3, Np = 150, seed = 11)
    ## scattered missingness in one continuous and one binary column
    set.seed(12)
    d$Y[sample(nrow(d$Y), 40), "c2"] <- NA
    d$Y[sample(nrow(d$Y), 40), "b2"] <- NA
    fit <- suppressWarnings(mvn_infer_mlm_mixed(d$Y, d$X, d$study,
      binary = d$binary, iter = 400, chains = 2, cores = 1, refresh = 0,
      seed = 3))
    fit_mixed <<- list(fit = fit, d = d)
  }
  fit_mixed
}

test_that("mvn_infer_mlm_mixed fits, hides per-record parameters", {
  skip_on_cran()
  f <- get_fit()
  fit <- f$fit
  ## per-record u / ymiss and per-study Omega_s are not saved by default
  expect_false(any(grepl("^(u|ymiss|Omega_s)\\[", names(fit))))
  expect_equal(fit@par_dims$Omega_global, c(4, 4))
  meta <- attr(fit, "muvamere")
  expect_equal(meta$var_names, c("c1", "b1", "c2", "b2"))
  expect_equal(meta$binary, f$d$binary)
  expect_s3_class(attr(fit, "diagnose"), "muvamere_diagnostics")
})

test_that("mixed-fit hyperparameters come back in the user's column order", {
  skip_on_cran()
  f <- get_fit()
  h <- mvn_extract_hyperparams(f$fit)
  expect_equal(colnames(h$OmegaG), c("c1", "b1", "c2", "b2"))
  expect_equal(colnames(h$betag), c("c1", "b1", "c2", "b2"))
  expect_equal(names(h$ltaum), c("c1", "c2")) # continuous variates only
  expect_equal(h$binary, c(FALSE, TRUE, FALSE, TRUE))
  ## the known c1-b1 latent (biserial) correlation of 0.6 is recovered, in
  ## the right cell despite the internal reordering
  expect_true(abs(h$OmegaG["c1", "b1"] - 0.6) < 0.2)
  expect_true(all(abs(h$OmegaG[upper.tri(h$OmegaG)][-1]) < 0.25))
  expect_true(min(eigen(h$OmegaG, TRUE, only.values = TRUE)$values) > 0)
})

test_that("mvn_infer_mlm_mixed runs with binary variates only", {
  skip_on_cran()
  set.seed(21)
  X <- cbind(1, runif(120))
  Y <- mvn_sample_study_kappa(X, matrix(0, 2, 3), matrix(0.05, 2, 3), 0.05,
    numeric(0), numeric(0), diag(3), binary = rep(TRUE, 3))
  Y <- rbind(Y, mvn_sample_study_kappa(X, matrix(0, 2, 3),
    matrix(0.05, 2, 3), 0.05, numeric(0), numeric(0), diag(3),
    binary = rep(TRUE, 3)))
  fit <- suppressWarnings(mvn_infer_mlm_mixed(Y, rbind(X, X),
    rep(1:2, each = 120), binary = 1:3, iter = 100, chains = 1, cores = 1,
    refresh = 0, seed = 1))
  expect_equal(fit@par_dims$ltaum, 0)
  h <- mvn_extract_hyperparams(fit)
  expect_length(h$ltaum, 0)
  expect_equal(dim(h$OmegaG), c(3, 3))
})

test_that("a generator holds no data, round-trips, and simulates", {
  skip_on_cran()
  f <- get_fit()
  g <- mvn_make_generator(f$fit, ndraws = 100)
  expect_s3_class(g, "muvamere_generator")
  expect_equal(g$ndraws, 100)
  expect_true(mvn_check_generator(g, n_records = nrow(f$d$Y)))
  expect_equal(dim(g$draws$Omega_global), c(100, 4, 4))
  expect_equal(dimnames(g$draws$Omega_global)[[2]], c("c1", "b1", "c2", "b2"))
  ## no stanfit, environment or record-level array anywhere
  sz <- as.numeric(object.size(g))
  expect_lt(sz, as.numeric(object.size(f$fit)) / 5)
  tmp <- tempfile(fileext = ".rds")
  saveRDS(g, tmp)
  g2 <- readRDS(tmp)
  expect_equal(g2, g)
  set.seed(8)
  AP <- mvn_generate_AP(g2, list(cbind(1, runif(50)), cbind(1, runif(30))))
  expect_equal(names(AP), c("c1", "b1", "c2", "b2", "studyno", "obsno"))
  expect_equal(nrow(AP), 80)
  expect_true(all(unlist(AP[, c("b1", "b2")]) %in% c(0, 1)))
  expect_false(all(AP$c1 %in% c(0, 1)))
  ## a fit works too (posterior-mean plug-in)
  AP2 <- mvn_generate_AP(f$fit, list(cbind(1, runif(10))))
  expect_equal(names(AP2)[1:4], c("c1", "b1", "c2", "b2"))
  expect_error(mvn_generate_AP(g, list(cbind(1, runif(5), 0))),
    "covariates")
})

test_that("mvn_check_generator catches data and foreign objects", {
  skip_on_cran()
  g <- mvn_make_generator(get_fit()$fit, ndraws = 50)
  n <- nrow(get_fit()$d$Y)
  bad <- g; bad$records <- 1:3
  expect_error(mvn_check_generator(bad), "unexpected components")
  bad <- g; bad$var_names <- c(bad$var_names, rep("x", n - 4))
  expect_error(mvn_check_generator(bad, n_records = n), "n_records")
  bad <- g; bad$draws$u <- matrix(0, 50, n)
  expect_error(mvn_check_generator(bad), "unexpected draws")
  bad <- g; bad$diagnostics$fit <- get_fit()$fit
  expect_error(mvn_check_generator(bad), "not allowed")
  bad <- g; bad$draws$Omega_global <- bad$draws$Omega_global[1:10, , ]
  expect_error(mvn_check_generator(bad), "has dim")
})

test_that("a generator can be made from an all-continuous fit", {
  skip_on_cran()
  set.seed(31)
  X <- cbind(1, runif(100))
  Ys <- lapply(1:2, function(i) mvn_sample_study_kappa(X,
    matrix(1, 2, 3), matrix(0.05, 2, 3), 0.05, rep(0, 3), rep(0.1, 3),
    diag(3)))
  fit <- suppressWarnings(mvn_infer_mlm_sparse(do.call(rbind, Ys),
    rbind(X, X), rep(1:2, each = 100), iter = 200, chains = 1, cores = 1,
    refresh = 0, seed = 2))
  g <- mvn_make_generator(fit, ndraws = 20)
  expect_equal(g$binary, rep(FALSE, 3))
  expect_true(mvn_check_generator(g, n_records = 200))
  AP <- mvn_generate_AP(g, list(cbind(1, runif(15))))
  expect_equal(dim(AP), c(15, 5))
})
