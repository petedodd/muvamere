## tests for ordinal variates in the mixed model: simulator, input
## handling, fitting, and the generator round trip.

## mixed data with an ordinal variate: user order (c1, o1, b1), o1 with 4
## levels (cutpoints 0, 0.6, 1.2 on its latent scale), c1-o1 latent
## (polyserial) correlation r
sim_ord <- function(S = 3, Np = 200, seed = 1, r = 0.6) {
  set.seed(seed)
  OG <- diag(3)
  OG[1, 2] <- OG[2, 1] <- r
  betag <- matrix(c(1, 0.2, 0.3, 0.2, -0.2, 0.3), 2, 3)
  colnames(betag) <- c("c1", "o1", "b1")
  Xlist <- replicate(S, cbind(1, runif(Np)), simplify = FALSE)
  Ys <- lapply(Xlist, function(X) {
    mvn_sample_study_kappa(X, betag, matrix(0.05, 2, 3),
      kappa = 0.05, ltaum = 0, lsig = 0.1, OmegaG = OG,
      binary = c(FALSE, FALSE, TRUE), cuts = list(o1 = c(0.6, 1.2))
    )
  })
  Y <- as.data.frame(do.call(rbind, Ys))
  Y$o1 <- factor(c("none", "mild", "moderate", "severe")[Y$o1 + 1],
    levels = c("none", "mild", "moderate", "severe"), ordered = TRUE
  )
  list(
    Y = Y, X = do.call(rbind, Xlist),
    study = rep(seq_len(S), each = Np), OG = OG
  )
}

test_that("mvn_sample_study_kappa emits ordinal codes on the latent scale", {
  set.seed(4)
  X <- matrix(1, 20000, 1)
  betag <- matrix(0, 1, 2, dimnames = list(NULL, c("a", "o")))
  Y <- mvn_sample_study_kappa(X, betag, matrix(0, 1, 2), kappa = 0,
    ltaum = 0, lsig = 0, OmegaG = diag(2), cuts = list(o = 0.5))
  expect_true(all(Y[, "o"] %in% 0:2))
  p <- tabulate(Y[, "o"] + 1, 3) / nrow(Y)
  expect_equal(p, c(0.5, pnorm(0.5) - 0.5, 1 - pnorm(0.5)),
    tolerance = 0.02)
  expect_error(mvn_sample_study_kappa(X, betag, matrix(0, 1, 2), 0, 0, 0,
    diag(2), cuts = list(z = 0.5)), "named by columns")
  expect_error(mvn_sample_study_kappa(X, betag, matrix(0, 1, 2), 0,
    numeric(0), numeric(0), diag(2), binary = c(FALSE, TRUE),
    cuts = list(o = 0.5)), "binary and ordinal")
  expect_error(mvn_sample_study_kappa(X, betag, matrix(0, 1, 2), 0, 0, 0,
    diag(2), cuts = list(o = -1)), "increasing")
})

test_that("ordinal input is encoded and checked before sampling", {
  enc <- muvamere:::.encode_ordinal(
    data.frame(a = c(3, 1, 2, NA), b = factor(c("x", "y", "z", "x"),
      levels = c("z", "y", "x"), ordered = TRUE)),
    c(TRUE, TRUE)
  )
  expect_equal(enc$Y$a, c(2L, 0L, 1L, NA))
  expect_equal(enc$Y$b, c(2L, 1L, 0L, 2L))
  expect_equal(enc$levels, list(a = c(1, 2, 3), b = c("z", "y", "x")))
  d <- sim_ord(S = 2, Np = 60)
  expect_error(mvn_infer_mlm_mixed(d$Y, d$X, d$study, binary = "b1",
    ordinal = c("o1", "b1")), "both binary and ordinal")
  expect_error(mvn_infer_mlm_mixed(d$Y, d$X, d$study, ordinal = "b1"),
    "use 'binary'")
  Y <- d$Y
  Y$o1[Y$o1 == "severe"] <- "moderate"
  Y$o1[1:2] <- "severe"
  expect_error(mvn_infer_mlm_mixed(Y, d$X, d$study, binary = "b1",
    ordinal = "o1"), "min_level")
})

## one shared small fit (fitting is the slow part)
fit_ord <- NULL
get_ord_fit <- function() {
  if (is.null(fit_ord)) {
    d <- sim_ord(S = 3, Np = 200, seed = 11)
    set.seed(12)
    d$Y$o1[sample(nrow(d$Y), 40)] <- NA # MAR (here MCAR) ordinal cells
    fit <- suppressWarnings(mvn_infer_mlm_mixed(d$Y, d$X, d$study,
      binary = "b1", ordinal = "o1", iter = 400, chains = 2, cores = 1,
      refresh = 0, seed = 3))
    fit_ord <<- list(fit = fit, d = d)
  }
  fit_ord
}

test_that("an ordinal fit recovers cutpoints and the latent correlation", {
  skip_on_cran()
  f <- get_ord_fit()
  meta <- attr(f$fit, "muvamere")
  expect_equal(meta$nlev, c(0L, 4L, 2L))
  expect_equal(meta$levels$o1, c("none", "mild", "moderate", "severe"))
  expect_equal(f$fit@par_dims$cuts, 2)
  h <- mvn_extract_hyperparams(f$fit)
  expect_equal(h$cuts$o1, c(0.6, 1.2), tolerance = 0.2)
  expect_true(abs(h$OmegaG["c1", "o1"] - 0.6) < 0.2)
  expect_true(abs(h$OmegaG["c1", "b1"]) < 0.25)
  expect_length(h$ltaum, 1) # c1 only
})

test_that("an ordinal generator checks, round-trips and returns levels", {
  skip_on_cran()
  f <- get_ord_fit()
  g <- mvn_make_generator(f$fit, ndraws = 100)
  expect_equal(dim(g$draws$cuts$o1), c(100, 2))
  expect_true(mvn_check_generator(g, n_records = nrow(f$d$Y)))
  tmp <- tempfile(fileext = ".rds")
  saveRDS(g, tmp)
  g2 <- readRDS(tmp)
  expect_equal(g2, g)
  set.seed(9)
  ## many synthetic new studies: each draws its own coefficients, so one
  ## study alone can differ from the pooled data by more than sampling noise
  AP <- mvn_generate_AP(g2, replicate(40, cbind(1, runif(100)),
    simplify = FALSE))
  expect_true(is.ordered(AP$o1))
  expect_equal(levels(AP$o1), c("none", "mild", "moderate", "severe"))
  expect_true(all(AP$b1 %in% c(0, 1)))
  ## pooled level frequencies of synthetic vs fitted data roughly agree
  p_dat <- prop.table(table(f$d$Y$o1))
  p_syn <- prop.table(table(AP$o1))
  expect_lt(max(abs(p_dat - p_syn)), 0.07)
  AP2 <- mvn_generate_AP(f$fit, list(cbind(1, runif(10))))
  expect_true(is.ordered(AP2$o1))
  bad <- g
  bad$draws$cuts$o1 <- bad$draws$cuts$o1[, 1, drop = FALSE]
  expect_error(mvn_check_generator(bad), "cuts\\$o1 has dim")
  bad <- g
  bad$draws$cuts <- NULL
  expect_error(mvn_check_generator(bad), "exactly the ordinal")
})

test_that("generators made before ordinal support still work", {
  skip_on_cran()
  set.seed(31)
  X <- cbind(1, runif(100))
  Ys <- lapply(1:2, function(i) mvn_sample_study_kappa(X,
    matrix(0, 2, 3), matrix(0.05, 2, 3), 0.05, rep(0, 2), rep(0.1, 2),
    diag(3), binary = c(FALSE, FALSE, TRUE)))
  fit <- suppressWarnings(mvn_infer_mlm_mixed(do.call(rbind, Ys),
    rbind(X, X), rep(1:2, each = 100), binary = 3, iter = 100,
    chains = 1, cores = 1, refresh = 0, seed = 2))
  g <- mvn_make_generator(fit, ndraws = 20)
  old <- g
  old$nlev <- old$levels <- NULL
  old$draws$cuts <- NULL
  expect_true(mvn_check_generator(old))
  AP <- mvn_generate_AP(old, list(cbind(1, runif(15))))
  expect_equal(dim(AP), c(15, 5))
  expect_true(all(AP$V3 %in% c(0, 1)))
})

test_that("the ordinal GHK step is the u-quantile, continuous in a", {
  skip_on_cran()
  ## compile only the functions block of the installed mixed model
  code <- stanmodels$mvn_inferH_mixed@model_code
  fn <- sub("(?s)^.*?(functions\\s*\\{.*?\n\\}).*$", "\\1", code,
    perl = TRUE)
  e <- new.env()
  rstan::expose_stan_functions(rstan::stanc(model_code = fn), env = e)
  ghk <- e$ghk_interval
  ## z must be the u-quantile of N(0,1) truncated to (a, b), and
  ## log P = log(Phi(b) - Phi(a)), on both sides of the a = 0 branch
  ## switch (a regression test: the a > 0 branch once drew the
  ## (1 - u)-quantile) and in the tails
  for (ab in list(c(-1e-9, 1), c(1e-9, 1), c(-2, -1), c(0.5, 0.8),
    c(3, 4), c(-1, 2))) {
    for (u in c(0.05, 0.5, 0.95)) {
      ## lo = a, hi = b with pre = 0 and Ljj = 1
      g <- ghk(0, 1, 1, ab[1], 1, ab[2], u)
      p <- pnorm(ab[2]) - pnorm(ab[1])
      expect_equal(g[1], log(p), tolerance = 1e-8)
      expect_equal(g[2], qnorm(pnorm(ab[1]) + p * u), tolerance = 1e-6)
    }
  }
  ## one-sided intervals match the binary code: bottom = u-quantile, top
  ## = (1 - u)-quantile, as in the original GHK steps
  expect_equal(ghk(0, 1, 0, 0, 1, -0.3, 0.2)[2], qnorm(pnorm(-0.3) * 0.2))
  expect_equal(ghk(0, 1, 1, 0.4, 0, 0, 0.2)[2],
    -qnorm(pnorm(-0.4) * 0.2))
})
