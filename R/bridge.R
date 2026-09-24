## this file bridges Stan inference (mvn_infer_mlm_sparse, mvn_infer_mlm_mixed)
## to the R-side simulator (mvn_sample_study_kappa), which otherwise expects
## its hyperparameter arguments to be supplied by hand.


##' Extract global hyperparameter point estimates from a fitted hierarchical
##' model
##'
##' Pulls posterior means (not full posterior draws) of the global
##' hyperparameters from a stanfit returned by \code{mvn_infer_mlm_sparse()}
##' or \code{mvn_infer_mlm_mixed()}, in the original variable order and
##' named, matching the arguments \code{mvn_sample_study_kappa()} expects.
##' This is a point-estimate plug-in: it collapses the posterior to its mean,
##' so a synthetic population built from these estimates will understate
##' parameter uncertainty; \code{mvn_make_generator()} keeps posterior draws
##' instead.
##'
##' @title mvn_extract_hyperparams
##' @param fit a stanfit from \code{mvn_infer_mlm_sparse()} or
##'   \code{mvn_infer_mlm_mixed()}
##' @return a list with elements \code{betag} (NP x NV matrix, posterior mean
##'   of BetaM), \code{sigb} (NP x NV matrix, posterior mean of BetaS),
##'   \code{ltaum}/\code{lsig} (posterior means of the mean and SD of
##'   log(tau) across studies, one per \emph{continuous} variate),
##'   \code{OmegaG} (NV x NV matrix, posterior mean of Omega_global),
##'   \code{model} (\code{"kappa"}), \code{kappa} (scalar) and \code{binary}
##'   (logical, which variates are binary)
##' @author Pete Dodd
##' @export
mvn_extract_hyperparams <- function(fit) {
  d <- .global_draws(fit)
  meta <- .meta_of(fit)
  m3 <- function(a) {
    out <- apply(a, c(2, 3), mean)
    dimnames(out) <- dimnames(a)[2:3]
    out
  }
  list(
    betag = m3(d$BetaM),
    sigb = m3(d$BetaS),
    ltaum = colMeans(d$ltaum),
    lsig = colMeans(d$lsig),
    OmegaG = m3(d$Omega_global),
    model = "kappa",
    kappa = mean(d$kappa),
    binary = meta$binary
  )
}


##' Generate a synthetic multi-study population from a fitted hierarchical
##' model
##'
##' Simulates one new study per target design matrix in \code{Xlist}, each
##' with its own study-level correlations, regression coefficients and
##' scales drawn around the fitted global hyperparameters (so the generated
##' cohorts reflect the fitted between-study heterogeneity), via
##' \code{mvn_sample_study_kappa()}. Binary variates are returned as 0/1.
##'
##' Given a \emph{generator} (\code{mvn_make_generator()}), each synthetic
##' study uses a different random posterior draw of the global
##' hyperparameters, so posterior uncertainty is propagated. Given a
##' \emph{fit}, the posterior means are plugged in for every study
##' (\code{mvn_extract_hyperparams()}), which understates that uncertainty.
##'
##' @title mvn_generate_AP
##' @param fit a \code{muvamere_generator}, or a stanfit from
##'   \code{mvn_infer_mlm_sparse()} or \code{mvn_infer_mlm_mixed()}
##' @param Xlist list of covariate data matrices, one per synthetic study to
##'   generate, with the same columns as the covariates used in the fit
##' @return a data frame with responses for all generated studies, columns
##'   named as the fitted variates, with studyno and obsno columns appended
##'   (same format as \code{mvn_simulate_studies()})
##' @author Pete Dodd
##' @export
mvn_generate_AP <- function(fit, Xlist) {
  gen <- inherits(fit, "muvamere_generator")
  if (gen) {
    d <- fit$draws
    binary <- fit$binary
    nms <- fit$var_names
    np <- length(fit$x_names)
  } else {
    hyper <- mvn_extract_hyperparams(fit)
    binary <- hyper$binary
    nms <- colnames(hyper$betag)
    np <- nrow(hyper$betag)
  }
  SS <- list()
  for (i in seq_along(Xlist)) {
    if (ncol(Xlist[[i]]) != np) {
      stop("Xlist[[", i, "]] has ", ncol(Xlist[[i]]),
        " columns; the fit has ", np, " covariates")
    }
    if (gen) {
      k <- sample.int(fit$ndraws, 1)
      betag <- d$BetaM[k, , , drop = TRUE]
      sigb <- d$BetaS[k, , , drop = TRUE]
      if (np == 1) {
        betag <- matrix(betag, 1)
        sigb <- matrix(sigb, 1)
      }
      colnames(betag) <- nms
      SS[[i]] <- mvn_sample_study_kappa(
        Xlist[[i]], betag, sigb, d$kappa[k], d$ltaum[k, ], d$lsig[k, ],
        d$Omega_global[k, , ], binary = binary
      )
    } else {
      SS[[i]] <- mvn_sample_study_kappa(
        Xlist[[i]], hyper$betag, hyper$sigb,
        hyper$kappa, hyper$ltaum, hyper$lsig, hyper$OmegaG,
        binary = binary
      )
    }
    SS[[i]] <- as.data.frame(SS[[i]])
    SS[[i]]$studyno <- i
    SS[[i]]$obsno <- 1:nrow(SS[[i]])
  }
  do.call("rbind", SS)
}
