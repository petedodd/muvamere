## the mixed continuous + binary hierarchical model (multivariate probit
## observation layer on the deviation-penalty model), and the metadata that
## lets downstream functions map the model's internal variable order back to
## the user's columns.

## Metadata attached to every hierarchical fit as attr(fit, "muvamere"):
## variable names and types in the USER's column order, the permutation to
## the model's internal order (continuous first), covariate names, and which
## Stan model was used. No record-level information.
.fit_meta <- function(Y, X, binary, model) {
  perm <- c(which(!binary), which(binary))
  list(
    model = model,
    var_names = colnames(Y),
    binary = unname(binary),
    perm = perm, # internal variable k = user column perm[k]
    x_names = colnames(X)
  )
}

## default column names where the user gave none
.named <- function(M, prefix) {
  M <- as.matrix(M)
  if (is.null(colnames(M))) colnames(M) <- paste0(prefix, seq_len(ncol(M)))
  M
}

## which columns are binary: logical, integer indices or column names
.binary_cols <- function(binary, Y) {
  NV <- ncol(Y)
  if (is.logical(binary)) {
    if (length(binary) != NV) {
      stop("logical 'binary' must have ncol(Y) entries")
    }
    return(binary)
  }
  if (is.character(binary)) {
    bad <- setdiff(binary, colnames(Y))
    if (length(bad)) {
      stop("'binary' names not in Y: ", paste(bad, collapse = ", "))
    }
    return(colnames(Y) %in% binary)
  }
  if (is.numeric(binary)) {
    if (any(binary < 1 | binary > NV)) stop("'binary' indices out of range")
    return(seq_len(NV) %in% binary)
  }
  stop("'binary' must be logical, column indices or column names")
}

## Guard rail: every variable, and every pair of variables, must be jointly
## observed in at least min_joint records within EVERY study. A variable (or
## pair) with no data in some study leaves that study's correlations
## unidentified; in testing this was the main source of divergences, and
## supporting it (per-study marginalisation) is out of scope.
.check_joint_observation <- function(obs, study, min_joint, var_names) {
  problems <- character()
  for (s in sort(unique(study))) {
    o <- obs[study == s, , drop = FALSE]
    joint <- crossprod(o * 1) # [NV, NV]: diagonal = per-variable counts
    low <- which(joint < min_joint & upper.tri(joint, diag = TRUE),
      arr.ind = TRUE)
    for (k in seq_len(nrow(low))) {
      i <- low[k, 1]; j <- low[k, 2]
      what <- if (i == j) var_names[i] else
        paste0(var_names[i], " & ", var_names[j])
      problems <- c(problems, sprintf("study %s: %s (%d jointly observed)",
        s, what, joint[i, j]))
    }
  }
  if (length(problems)) {
    stop("too few jointly observed records (min_joint = ", min_joint, "):\n",
      paste(" ", utils::head(problems, 20), collapse = "\n"),
      if (length(problems) > 20) "\n  ..." else "",
      "\nVariables missing (almost) entirely from a study are not supported;",
      " drop the variable or the study.", call. = FALSE)
  }
  invisible(TRUE)
}


##' MCMC sampling for mixed continuous and binary data from multiple studies
##'
##' The mixed-type version of \code{mvn_infer_mlm_sparse()}: the same
##' hierarchical deviation-penalty model (per-study regressions partially
##' pooled; a regularized-horseshoe global correlation matrix; each study's
##' correlations deviating from it with estimated strength \code{kappa}),
##' with binary variates modelled by a multivariate \emph{probit} layer. Each
##' binary is the indicator that a latent normal variate exceeds 0; the
##' latent variates have unit scale, so the correlation matrix holds latent
##' (tetrachoric, for two binaries; biserial, for a binary and a continuous
##' variate) correlations, which are fully identified.
##'
##' \strong{Missing data.} Missing values (\code{NA}) are allowed in
##' \code{Y}, in both continuous and binary columns, and are handled exactly
##' under missing-at-random: a missing binary drops out of that record's
##' likelihood; a missing continuous value is estimated as a parameter.
##' \code{X} must be complete. Every variable, and every pair of variables,
##' must be jointly observed in at least \code{min_joint} records in every
##' study: a variable missing (almost) entirely from a study is not supported
##' and gives an error.
##'
##' \strong{Privacy.} The model has per-record parameters (one GHK uniform
##' per record and binary variate, and the imputed missing continuous
##' values) that, together with the other parameters, reveal individual
##' records. They are not saved in the returned fit unless \code{pars} is
##' given; for sharing, use \code{mvn_make_generator()}.
##'
##' \strong{Cost.} The likelihood is evaluated record by record for the
##' binary variates, so fits are slower than the all-continuous model: about
##' an hour for 2 chains x 1000 iterations with 10 variates and ~2200 records
##' on a laptop in testing, and much longer for 20 variates.
##'
##' @title mvn_infer_mlm_mixed
##' @param Y outcome data (matrix or data frame), one row per record: numeric
##'   continuous columns and binary columns coded 0/1 (or logical), with
##'   \code{NA} for missing values
##' @param X covariate matrix, one row per record, no missing values
##' @param study vector of which study each record belongs to
##' @param binary which columns of \code{Y} are binary: a logical vector,
##'   column indices, or column names. At least one is required (use
##'   \code{mvn_infer_mlm_sparse()} for all-continuous data).
##' @param betaM_prior_sd prior SD for the global regression coefficient means
##' @param betaS_prior_sd prior SD (half-normal) for the between-study
##'   regression coefficient SDs
##' @param tauM_prior_sd prior SD for the mean of log(tau) across studies
##'   (continuous variates; binary latent scales are fixed at 1)
##' @param tauS_prior_sd prior SD (half-normal) for the SD of log(tau) across
##'   studies
##' @param p0 prior guess for the number of non-zero correlations, as in
##'   \code{mvn_infer_mlm_sparse()}
##' @param slab_scale slab scale for large correlations, default 0.5
##' @param slab_df slab degrees of freedom, default 4
##' @param kappa_prior_scale scale of the half-normal prior on \code{kappa},
##'   default 2
##' @param min_joint minimum number of records in which every variable and
##'   every pair of variables must be jointly observed, in every study
##'   (default 5)
##' @param init optional Stan \code{init}; default \code{NULL} gives informed
##'   starts
##' @param iter iterations for MCMC
##' @param cores number of cores to use
##' @param chains number of chains to use
##' @param ... further arguments passed to \code{rstan::sampling()}. Unless
##'   \code{pars} is given, the per-record \code{u} and \code{ymiss} and the
##'   per-study \code{Omega_s} are not saved.
##' @return a Stan sample object with attributes \code{"diagnose"} (the
##'   \code{mvn_diagnose()} result) and \code{"muvamere"} (variable names and
##'   types, covariate names; no record-level data). Parameters are in the
##'   model's internal order (continuous variates first);
##'   \code{mvn_extract_hyperparams()} and \code{mvn_make_generator()} map
##'   them back to the columns of \code{Y}.
##' @author Pete Dodd
##' @export
##' @import rstan
mvn_infer_mlm_mixed <- function(Y, X, study, binary,
                                betaM_prior_sd = 1,
                                betaS_prior_sd = 0.5,
                                tauM_prior_sd = 1,
                                tauS_prior_sd = 0.5,
                                p0 = NULL, slab_scale = 0.5, slab_df = 4,
                                kappa_prior_scale = 2,
                                min_joint = 5,
                                init = NULL,
                                iter = 2e3, cores = 4, chains = 4, ...) {
  Y <- .named(Y, "V")
  X <- .named(X, "X")
  binary <- .binary_cols(binary, Y)
  if (!any(binary)) {
    stop("no binary columns: use mvn_infer_mlm_sparse() for continuous data")
  }
  if (ncol(Y) < 2) stop("need at least 2 variates (columns of Y)")
  if (anyNA(X)) stop("X must not contain missing values")
  Yb_user <- Y[, binary, drop = FALSE]
  if (is.logical(Yb_user)) Yb_user <- Yb_user * 1
  if (any(!is.na(Yb_user) & !(Yb_user %in% c(0, 1)))) {
    stop("binary columns must contain only 0, 1 (or TRUE/FALSE) and NA")
  }
  meta <- .fit_meta(Y, X, binary, "mixed")
  ## internal order: continuous first, then binary
  Y <- Y[, meta$perm, drop = FALSE]
  NC <- sum(!binary)
  NB <- sum(binary)
  NV <- NC + NB
  srt <- .sort_by_study(Y, X, study)
  Y <- srt$Y
  X <- srt$X
  study <- srt$study
  obs <- !is.na(Y)
  .check_joint_observation(obs, study, min_joint, colnames(Y))
  ## a binary constant within a study carries no information on that
  ## study's correlations with it: allowed, but flagged
  for (j in NC + seq_len(NB)) {
    for (s in unique(study)) {
      v <- Y[study == s & obs[, j], j]
      if (length(unique(v)) == 1) {
        warning("binary '", colnames(Y)[j], "' is constant (", v[1],
          ") within study ", s, call. = FALSE)
      }
    }
  }
  DR <- choose(NV, 2)
  if (is.null(p0)) p0 <- .default_p0(NV)
  if (!(is.numeric(p0) && length(p0) == 1 && p0 > 0 && p0 < DR)) {
    stop("p0 must be a single number with 0 < p0 < choose(ncol(Y), 2) = ",
      DR)
  }
  if (!(slab_scale > 0 && slab_df > 0)) {
    stop("slab_scale and slab_df must be positive")
  }
  if (!(is.numeric(kappa_prior_scale) && length(kappa_prior_scale) == 1 &&
    kappa_prior_scale > 0)) {
    stop("kappa_prior_scale must be a single positive number")
  }
  Yc <- Y[, seq_len(NC), drop = FALSE]
  miss_c <- which(is.na(Yc), arr.ind = TRUE)
  miss_c <- matrix(as.integer(miss_c), ncol = 2)
  Yc[is.na(Yc)] <- 0 # any value: replaced by the ymiss parameters
  Yb <- Y[, NC + seq_len(NB), drop = FALSE]
  Yb[is.na(Yb)] <- -1
  storage.mode(Yb) <- "integer"
  shdata <- list(
    Nrecords = nrow(Y), Nstudies = length(unique(study)), study = study,
    NP = ncol(X), NC = NC, NB = NB, X = X, Yc = Yc, Yb = Yb,
    Nmiss_c = nrow(miss_c), miss_c = miss_c,
    betaM_prior_sd = betaM_prior_sd, betaS_prior_sd = betaS_prior_sd,
    tauM_prior_sd = tauM_prior_sd, tauS_prior_sd = tauS_prior_sd,
    p0 = p0, slab_scale = slab_scale, slab_df = slab_df,
    kappa_prior_scale = kappa_prior_scale
  )
  if (is.null(init)) {
    ## informed Omega_global start from the crude data (missing -> column
    ## mean); u = 0.5; missing continuous values start at column means
    crude <- apply(Y, 2, function(v) {
      v[is.na(v)] <- mean(v, na.rm = TRUE)
      v
    })
    cm <- colMeans(Y[, seq_len(NC), drop = FALSE], na.rm = TRUE)
    init <- lapply(.rhs_kappa_inits(crude, X, study, chains, slab_scale),
      function(l) {
        l$u <- matrix(0.5, nrow(Y), NB)
        if (nrow(miss_c)) {
          l$ymiss <- array(cm[miss_c[, 2]], nrow(miss_c))
        }
        l
      })
  }
  args <- list(...)
  if (is.null(args$pars)) {
    ## per-record u/ymiss reveal records; per-study Omega_s is unused
    args$pars <- c("u", "ymiss", "Omega_s")
    args$include <- FALSE
  }
  fit <- do.call(rstan::sampling, c(list(stanmodels$mvn_inferH_mixed,
    data = shdata, chains = chains, cores = cores, iter = iter,
    init = init
  ), args))
  fit <- .attach_stability_check(fit, "")
  attr(fit, "muvamere") <- meta
  fit
}
