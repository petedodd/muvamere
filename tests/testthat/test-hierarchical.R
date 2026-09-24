## tests for the hierarchical model's bookkeeping (mvn_infer_mlm_sparse()) and
## the R/bridge.R functions built on top of it. (Rewritten 2026-09-24 when
## the rho-blend mvn_infer_mlm() was removed; same checks, current model.)

## small multi-study fixture from the package's own simulator
sim_studies <- function(sizes, NV = 3, seed = 301) {
  set.seed(seed)
  Xlist <- lapply(sizes, function(n) cbind(1, runif(n)))
  allsamp <- mvn_simulate_studies(
    Xlist,
    betag = matrix(1, 2, NV), sigb = matrix(0.05, 2, NV),
    kappa = 0.05, ltaum = rep(0, NV), lsig = rep(0.1, NV), lkj_global = 2
  )
  list(Y = as.matrix(allsamp[, 1:NV]), X = do.call(rbind, Xlist),
    study = allsamp$studyno)
}

test_that("Nstudies comes from distinct studies, not Nrecords", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  ## Nrecords = 3 * 25 = 75 -- must NOT show up as a per-study dimension
  d <- sim_studies(rep(25, 3))
  expect_equal(nrow(d$Y), 75) # sanity check on the fixture itself
  fit <- suppressWarnings(mvn_infer_mlm_sparse(d$Y, d$X, d$study,
    iter = 300, chains = 1, cores = 1, refresh = 0
  ))
  expect_equal(fit@par_dims$tau, c(3, 3))
  expect_equal(fit@par_dims$rr, c(3, 3))
  expect_equal(fit@par_dims$Betas, c(3, 2, 3))
})

test_that("mvn_extract_hyperparams gives correctly (re)shaped hyperparams", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  d <- sim_studies(rep(25, 3), seed = 302)
  fit <- suppressWarnings(mvn_infer_mlm_sparse(d$Y, d$X, d$study,
    iter = 300, chains = 1, cores = 1, refresh = 0
  ))
  hyper <- mvn_extract_hyperparams(fit)
  expect_named(
    hyper,
    c("betag", "sigb", "ltaum", "lsig", "OmegaG", "model", "kappa", "binary")
  )
  expect_equal(hyper$model, "kappa")
  expect_equal(dim(hyper$betag), c(2, 3))
  expect_equal(dim(hyper$sigb), c(2, 3))
  expect_length(hyper$ltaum, 3)
  expect_length(hyper$lsig, 3)
  expect_length(hyper$kappa, 1)
  expect_true(hyper$kappa >= 0)
  expect_equal(dim(hyper$OmegaG), c(3, 3))
  ## a correlation matrix
  expect_equal(unname(diag(hyper$OmegaG)), rep(1, 3), tolerance = 1e-6)
  expect_equal(hyper$binary, rep(FALSE, 3))
  ## outputs are named after the columns of Y (default names here)
  expect_equal(colnames(hyper$OmegaG), c("V1", "V2", "V3"))
})

test_that("mvn_extract_hyperparams refuses a non-kappa fit", {
  ## minimal stand-in with a par_dims slot, as a (removed) rho-blend fit had
  methods::setClass("fakefit", representation(par_dims = "list"))
  fake <- methods::new("fakefit", par_dims = list(rho = integer(0)))
  expect_error(mvn_extract_hyperparams(fake), "no kappa parameter")
})

test_that("mvn_generate_AP output matches mvn_simulate_studies format", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  d <- sim_studies(rep(25, 3), seed = 303)
  fit <- suppressWarnings(mvn_infer_mlm_sparse(d$Y, d$X, d$study,
    iter = 300, chains = 1, cores = 1, refresh = 0
  ))
  ## generate for a different set of target studies
  Xtarget <- replicate(2, cbind(1, runif(10)), simplify = FALSE)
  AP <- mvn_generate_AP(fit, Xtarget)

  expect_s3_class(AP, "data.frame")
  expect_equal(nrow(AP), 20)
  expect_true(all(c("studyno", "obsno") %in% names(AP)))
  expect_equal(sort(unique(AP$studyno)), 1:2)
  expect_equal(as.vector(table(AP$studyno)), c(10, 10))
})


## the hierarchical Stan model evaluates the likelihood study-by-study and
## needs records sorted by study; the R wrapper sorts via .sort_by_study().

test_that(".sort_by_study makes contiguous study blocks, relabels ids 1..S", {
  Y <- cbind(a = 1:6, b = 11:16)
  X <- cbind(1, c(6, 5, 4, 3, 2, 1))
  study <- c(20, 10, 20, 30, 10, 30) # unsorted, ids not 1..S
  s <- muvamere:::.sort_by_study(Y, X, study)
  expect_equal(s$study, c(1L, 1L, 2L, 2L, 3L, 3L))
  ## rows moved with their study
  expect_equal(s$Y[, "a"], c(2L, 5L, 1L, 3L, 4L, 6L))
  expect_equal(s$X[, 2], X[order(study), 2]) # X moved consistently with Y
  expect_error(
    muvamere:::.sort_by_study(Y, X, study[-1]), "same number of records"
  )
})

test_that("record order does not change the mvn_infer_mlm_sparse fit", {
  skip_on_cran()
  skip_if_not_installed("trialr")
  d <- sim_studies(c(20, 35), NV = 2, seed = 304)
  st <- d$study
  ## interleave the studies but keep each study's own record order, so that
  ## after sorting the data are bitwise identical
  perm <- order(
    ave(seq_along(st), st, FUN = seq_along) + runif(length(st), 0, 0.5)
  )
  expect_true(is.unsorted(st[perm])) # fixture really is reordered
  ## identical inits: the wrapper would draw its jittered inits from R's RNG
  set.seed(9)
  ini <- muvamere:::.rhs_kappa_inits(d$Y, d$X, st, chains = 1,
    slab_scale = 0.5
  )
  fit1 <- suppressWarnings(mvn_infer_mlm_sparse(d$Y, d$X, st,
    init = ini, iter = 200, chains = 1, cores = 1, refresh = 0, seed = 5
  ))
  fit2 <- suppressWarnings(mvn_infer_mlm_sparse(
    d$Y[perm, ], d$X[perm, ], st[perm],
    init = ini, iter = 200, chains = 1, cores = 1, refresh = 0, seed = 5
  ))
  ## identical data after sorting, same seed => identical draws
  ## (permuted=FALSE: the default permuted=TRUE shuffles draws with R's RNG,
  ## which would differ)
  expect_equal(
    rstan::extract(fit1, pars = "kappa", permuted = FALSE),
    rstan::extract(fit2, pars = "kappa", permuted = FALSE)
  )
})
