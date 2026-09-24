## this file bridges Stan inference (mvn_infer_mlm_sparse) to the R-side
## simulator (mvn_sample_study_kappa), which otherwise expects its
## hyperparameter arguments to be supplied by hand.


##' Extract global hyperparameter point estimates from a fitted hierarchical
##' model
##'
##' Pulls posterior means (not full posterior draws) of the global
##' hyperparameters from a stanfit returned by \code{mvn_infer_mlm_sparse()},
##' reshaped to match the argument names \code{mvn_sample_study_kappa()}
##' expects. This is a point-estimate plug-in: it collapses the posterior to
##' its mean before generating, so a synthetic population built from these
##' estimates will understate parameter uncertainty relative to proper
##' posterior-predictive generation (drawing a fresh set of hyperparameters
##' per posterior draw, not just plugging in one point estimate).
##'
##' @title mvn_extract_hyperparams
##' @param fit a stanfit object returned by \code{mvn_infer_mlm_sparse()}
##' @return a list with elements \code{betag} (NP x NV matrix, posterior mean
##'   of BetaM), \code{sigb} (NP x NV matrix, posterior mean of BetaS),
##'   \code{ltaum}/\code{lsig} (length-NV vectors, posterior means of the mean
##'   and SD of log(tau) across studies), \code{OmegaG} (NV x NV matrix,
##'   posterior mean of Omega_global), \code{model} (\code{"kappa"}, the
##'   deviation-penalty model) and \code{kappa} (scalar)
##' @author Pete Dodd
##' @export
mvn_extract_hyperparams <- function(fit) {
  if (!"kappa" %in% names(fit@par_dims)) {
    stop("fit is not from mvn_infer_mlm_sparse() (no kappa parameter); ",
      "the rho-blend models were removed from muvamere in 2026-09")
  }
  NP <- fit@par_dims$BetaM[1]
  NV <- fit@par_dims$BetaM[2]

  ## NOTE reshape order: rstan::summary() lists a matrix parameter's elements
  ## in row-major order (par[1,1],par[1,2],...,par[2,1],...) not column-major,
  ## so byrow=TRUE is required to avoid a transpose
  get_mat <- function(par, nr, nc) {
    matrix(rstan::summary(fit, pars = par)$summary[, "mean"],
      nrow = nr, ncol = nc, byrow = TRUE
    )
  }
  get_vec <- function(par) {
    unname(rstan::summary(fit, pars = par)$summary[, "mean"])
  }

  list(
    betag = get_mat("BetaM", NP, NV),
    sigb = get_mat("BetaS", NP, NV),
    ltaum = get_vec("ltaum"),
    lsig = get_vec("lsig"),
    OmegaG = get_mat("Omega_global", NV, NV),
    model = "kappa",
    kappa = get_vec("kappa")
  )
}


##' Generate a synthetic multi-study population from a fitted hierarchical
##' model
##'
##' Plugs posterior-mean hyperparameters (via \code{mvn_extract_hyperparams()})
##' into \code{mvn_sample_study_kappa()} for each target design matrix in
##' \code{Xlist}. Each call draws a new study-level correlation and per-study
##' Betas/taus, so this generates new synthetic cohorts consistent with the
##' fitted between-study heterogeneity. Uses posterior means for the global
##' hyperparameters; it does not propagate posterior uncertainty in those
##' hyperparameters into the generated population. A fuller version should
##' draw hyperparameters from the posterior per replicate rather than plugging
##' in one point estimate (see the package TODO list).
##'
##' @title mvn_generate_AP
##' @param fit a stanfit object returned by \code{mvn_infer_mlm_sparse()}
##' @param Xlist list of covariate data matrices, one per synthetic study to
##'   generate
##' @return a data frame with responses for all generated studies, with studyno
##'   and obsno columns appended (same format as \code{mvn_simulate_studies()})
##' @author Pete Dodd
##' @export
mvn_generate_AP <- function(fit, Xlist) {
  hyper <- mvn_extract_hyperparams(fit)
  SS <- list()
  for (i in seq_along(Xlist)) {
    SS[[i]] <- mvn_sample_study_kappa(
      Xlist[[i]], hyper$betag, hyper$sigb,
      hyper$kappa, hyper$ltaum, hyper$lsig, hyper$OmegaG
    )
    SS[[i]] <- as.data.frame(SS[[i]])
    SS[[i]]$studyno <- i
    SS[[i]]$obsno <- 1:nrow(SS[[i]])
  }
  do.call("rbind", SS)
}
