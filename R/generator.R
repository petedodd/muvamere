## A data-free, saveable object for simulating new populations from a fitted
## hierarchical model. A stanfit can carry record-level information (e.g.
## the mixed model's per-record GHK uniforms); a generator holds only
## posterior draws of the global hyperparameters plus variable metadata.

## the variable-order metadata of a fit; fits made before the metadata
## attribute existed are treated as all-continuous with default names
.meta_of <- function(fit) {
  meta <- attr(fit, "muvamere")
  if (!is.null(meta)) {
    return(meta)
  }
  NV <- fit@par_dims$BetaM[2]
  NP <- fit@par_dims$BetaM[1]
  list(
    model = "kappa", var_names = paste0("V", seq_len(NV)),
    binary = rep(FALSE, NV), perm = seq_len(NV),
    x_names = paste0("X", seq_len(NP))
  )
}

## posterior draws of the global hyperparameters, selected by name and put
## in the user's variable order. Returns arrays with draws first:
## BetaM/BetaS [nd, NP, NV], ltaum/lsig [nd, NC], Omega_global [nd, NV, NV],
## kappa [nd].
.global_draws <- function(fit, idx = NULL) {
  if (!"kappa" %in% names(fit@par_dims)) {
    stop(
      "fit is not from mvn_infer_mlm_sparse() or mvn_infer_mlm_mixed() ",
      "(no kappa parameter); the rho-blend models were removed in 2026-09"
    )
  }
  meta <- .meta_of(fit)
  NP <- fit@par_dims$BetaM[1]
  NV <- fit@par_dims$BetaM[2]
  NC <- sum(!meta$binary)
  mat <- function(par, nr, nc) {
    A <- as.matrix(fit, pars = par)
    if (!is.null(idx)) A <- A[idx, , drop = FALSE]
    out <- array(NA_real_, c(nrow(A), nr, nc))
    for (i in seq_len(nr)) {
      for (j in seq_len(nc)) {
        out[, i, j] <- A[, sprintf("%s[%d,%d]", par, i, j)]
      }
    }
    out
  }
  vec <- function(par, n) {
    A <- as.matrix(fit, pars = par)
    if (!is.null(idx)) A <- A[idx, , drop = FALSE]
    A[, sprintf("%s[%d]", par, seq_len(n)), drop = FALSE]
  }
  back <- order(meta$perm) # internal -> user column order
  BetaM <- mat("BetaM", NP, NV)[, , back, drop = FALSE]
  BetaS <- mat("BetaS", NP, NV)[, , back, drop = FALSE]
  Og <- mat("Omega_global", NV, NV)[, back, back, drop = FALSE]
  k <- as.matrix(fit, pars = "kappa")[, 1]
  if (!is.null(idx)) k <- k[idx]
  nd <- length(k)
  ltaum <- if (NC) vec("ltaum", NC) else matrix(numeric(0), nd, 0)
  lsig <- if (NC) vec("lsig", NC) else matrix(numeric(0), nd, 0)
  dimnames(BetaM) <- dimnames(BetaS) <- list(
    NULL, meta$x_names,
    meta$var_names
  )
  dimnames(Og) <- list(NULL, meta$var_names, meta$var_names)
  cn <- meta$var_names[!meta$binary]
  dimnames(ltaum) <- dimnames(lsig) <- list(NULL, cn)
  list(
    BetaM = BetaM, BetaS = BetaS, ltaum = ltaum, lsig = lsig,
    Omega_global = Og, kappa = unname(k)
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a


##' Make a data-free generator object from a fitted hierarchical model
##'
##' Builds a small object that can simulate new populations (via
##' \code{mvn_generate_AP()}) but contains \strong{no study data}: only
##' posterior draws of the global hyperparameters (\code{BetaM},
##' \code{BetaS}, \code{ltaum}, \code{lsig}, \code{Omega_global},
##' \code{kappa}) in the original variable order, plus metadata (variable
##' and covariate names, which variates are binary, fit diagnostics and the
##' package version). Nothing record-level (e.g. the mixed model's GHK
##' uniforms or imputed missing values) and nothing study-level is kept, and
##' no stanfit, environment or function is stored, so the object can be
##' saved with \code{saveRDS()} and shared where the data cannot.
##' \code{mvn_check_generator()} audits an object for this.
##'
##' This removes the data, not all information derived from it: posterior
##' draws of hyperparameters are aggregate summaries of the fitted studies.
##' It is not a formal privacy guarantee (e.g. differential privacy).
##'
##' @title mvn_make_generator
##' @param fit a stanfit from \code{mvn_infer_mlm_sparse()} or
##'   \code{mvn_infer_mlm_mixed()}
##' @param ndraws number of posterior draws to keep (evenly thinned from the
##'   fit's draws; all of them if fewer are available), default 1000
##' @return an object of class \code{muvamere_generator}
##' @author Pete Dodd
##' @export
mvn_make_generator <- function(fit, ndraws = 1000) {
  total <- nrow(as.matrix(fit, pars = "kappa"))
  idx <- if (ndraws >= total) {
    seq_len(total)
  } else {
    unique(round(seq(1, total, length.out = ndraws)))
  }
  meta <- .meta_of(fit)
  dg <- attr(fit, "diagnose") %||% mvn_diagnose(fit)
  g <- list(
    draws = .global_draws(fit, idx),
    var_names = meta$var_names,
    binary = meta$binary,
    x_names = meta$x_names,
    model = meta$model,
    ndraws = length(idx),
    diagnostics = list(
      rhat_max = dg$rhat_max, ess_min = dg$ess_min,
      n_divergent = dg$n_divergent, n_draws = dg$n_draws, ok = dg$ok
    ),
    muvamere_version = as.character(utils::packageVersion("muvamere")),
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  class(g) <- "muvamere_generator"
  mvn_check_generator(g)
  g
}


##' Audit a generator object for anything that is not needed to simulate
##'
##' Checks that \code{g} has exactly the expected components, holds no
##' functions, environments or S4 objects (such as a stanfit), and that its
##' posterior-draw arrays have consistent dimensions. If \code{n_records}
##' (the number of records the model was fitted to) is given, also checks
##' that no component has a dimension or length equal to it: a signature
##' of record-level data. Used automatically by \code{mvn_make_generator()};
##' run it again before sharing an object from elsewhere.
##'
##' @title mvn_check_generator
##' @param g an object to check
##' @param n_records optional number of records in the fitted data
##' @return \code{invisible(TRUE)}, or an error listing the problems
##' @author Pete Dodd
##' @export
mvn_check_generator <- function(g, n_records = NULL) {
  top <- c(
    "draws", "var_names", "binary", "x_names", "model", "ndraws",
    "diagnostics", "muvamere_version", "created"
  )
  dr <- c("BetaM", "BetaS", "ltaum", "lsig", "Omega_global", "kappa")
  problems <- character()
  add <- function(...) problems <<- c(problems, paste0(...))
  if (!inherits(g, "muvamere_generator")) add("not a muvamere_generator")
  if (!is.list(g)) stop("not a list: ", class(g)[1])
  extra <- setdiff(names(g), top)
  if (length(extra)) add("unexpected components: ", toString(extra))
  extra <- setdiff(names(g$draws), dr)
  if (length(extra)) add("unexpected draws: ", toString(extra))
  miss <- setdiff(dr, names(g$draws))
  if (length(miss)) add("missing draws: ", toString(miss))
  ## no code, environments or S4 objects anywhere
  walk <- function(x, path) {
    if (is.function(x) || is.environment(x) || isS4(x)) {
      add(path, ": ", class(x)[1], " not allowed")
    } else if (is.list(x)) {
      for (nm in names(x)) walk(x[[nm]], paste0(path, "$", nm))
    }
  }
  walk(unclass(g), "g")
  ## draw arrays: leading dimension = ndraws, others consistent
  nd <- g$ndraws
  NV <- length(g$var_names)
  NC <- sum(!g$binary)
  NP <- length(g$x_names)
  want <- list(
    BetaM = c(nd, NP, NV), BetaS = c(nd, NP, NV),
    ltaum = c(nd, NC), lsig = c(nd, NC), Omega_global = c(nd, NV, NV)
  )
  for (nm in intersect(names(want), names(g$draws))) {
    d <- dim(g$draws[[nm]])
    if (!identical(as.numeric(d), as.numeric(want[[nm]]))) {
      add(
        "draws$", nm, " has dim ", toString(d), ", expected ",
        toString(want[[nm]])
      )
    }
  }
  if (length(g$draws$kappa) != nd) add("draws$kappa length != ndraws")
  ## nothing record-sized (the draws dimension itself is exempt)
  if (!is.null(n_records)) {
    sizes <- function(x, lead) {
      d <- dim(x) %||% length(x)
      if (lead && length(d) > 1) d[-1] else if (lead) integer(0) else d
    }
    for (nm in names(g$draws)) {
      if (any(sizes(g$draws[[nm]], TRUE) == n_records)) {
        add("draws$", nm, " has a dimension equal to n_records")
      }
    }
    for (nm in setdiff(names(g), c("draws", "ndraws"))) {
      x <- g[[nm]]
      xs <- if (is.list(x)) x else list(x)
      for (el in xs) {
        if (any(sizes(el, FALSE) == n_records)) {
          add(nm, " has a dimension/length equal to n_records")
        }
      }
    }
  }
  if (length(problems)) {
    stop("generator check failed:\n", paste(" ", problems, collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

##' @export
print.muvamere_generator <- function(x, ...) {
  cat("muvamere generator (", x$model, " model): ", length(x$var_names),
    " variates (", sum(x$binary), " binary), ", length(x$x_names),
    " covariates, ", x$ndraws, " posterior draws\n",
    sep = ""
  )
  cat("  variates:  ", paste(x$var_names, collapse = ", "), "\n", sep = "")
  cat("  covariates:", paste(x$x_names, collapse = ", "), "\n")
  cat("  fit: Rhat max ", signif(x$diagnostics$rhat_max, 3), ", ESS min ",
    signif(x$diagnostics$ess_min, 3), ", ", x$diagnostics$n_divergent,
    " divergent; muvamere ", x$muvamere_version, "\n",
    sep = ""
  )
  invisible(x)
}
