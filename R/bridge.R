## this file bridges Stan inference (mvn_infer_mlm) to the R-side simulators
## (mvn_sample_study / mvn_simulate_studies), which otherwise expect their
## hyperparameter arguments to be supplied by hand.



##' Extract global hyperparameter point estimates from a fitted hierarchical model
##'
##' Pulls posterior means (not full posterior draws) of the global
##' hyperparameters from a stanfit returned by \code{mvn_infer_mlm()}, reshaped
##' to match the argument names \code{mvn_sample_study()} already expects.
##' This is a point-estimate plug-in: it collapses the posterior to its mean
##' before generating, so a synthetic population built from these estimates
##' will understate parameter uncertainty relative to proper posterior-predictive
##' generation (drawing a fresh set of hyperparameters per posterior draw,
##' not just plugging in one point estimate).
##'
##' @title mvn_extract_hyperparams
##' @param fit a stanfit object returned by \code{mvn_infer_mlm()}
##' @return a list with elements \code{betag} (NP x NV matrix, posterior mean of BetaM),
##'   \code{sigb} (NP x NV matrix, posterior mean of BetaS), \code{taug} (length-NV vector,
##'   posterior mean of taum), \code{sigt} (length-NV vector, posterior mean of sigt),
##'   \code{rho} (scalar, posterior mean of rho), \code{OmegaG} (NV x NV matrix,
##'   posterior mean of Omega_global)
##' @author Pete Dodd
##' @export
mvn_extract_hyperparams <- function(fit) {
  NP <- fit@par_dims$BetaM[1]
  NV <- fit@par_dims$BetaM[2]

  ## NOTE reshape order: rstan::summary() lists a matrix parameter's elements
  ## in ROW-major order (par[1,1],par[1,2],...,par[2,1],...) not column-major,
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
    taug = get_vec("taum"),
    sigt = get_vec("sigt"),
    rho = get_vec("rho"),
    OmegaG = get_mat("Omega_global", NV, NV)
  )
}


##' Generate a synthetic multi-study population from a fitted hierarchical model
##'
##' Plugs posterior-mean hyperparameters (via \code{mvn_extract_hyperparams()})
##' into \code{mvn_sample_study()} for each target design matrix in
##' \code{Xlist}. Each call draws a new study-level local correlation and
##' per-study Betas/taus. So this generates new synthetic cohorts
##' consistent with the fitted between-study heterogeneity.
##' Uses posterior means for the global
##' hyperparameters (rho, Omega_global, BetaM, BetaS, taum, sigt); it does
##' not propagate posterior uncertainty in those hyperparameters into the
##' generated population. A fuller version should draw hyperparameters
##' from the posterior per replicate rather than plugging in one point
##' estimate (see the package TODO list).
##'
##' @title mvn_generate_AP
##' @param fit a stanfit object returned by \code{mvn_infer_mlm()}
##' @param Xlist list of covariate data matrices, one per synthetic study to generate
##' @param lkj_local local correlation LKJ prior parameter for the new studies'
##'   local correlation draws. This is a prior tuning constant that
##'   \code{mvn_infer_mlm()} takes as fixed input data (not something it
##'   estimates a posterior for) -- reuse the value passed to the original fit.
##' @return a data frame with responses for all generated studies, with studyno
##'   and obsno columns appended (same format as \code{mvn_simulate_studies()})
##' @author Pete Dodd
##' @export
mvn_generate_AP <- function(fit, Xlist, lkj_local = 3) {
  hyper <- mvn_extract_hyperparams(fit)
  SS <- list()
  for (i in seq_along(Xlist)) {
    SS[[i]] <- mvn_sample_study(
      Xlist[[i]], hyper$betag, hyper$sigb,
      hyper$rho, hyper$taug, hyper$sigt,
      hyper$OmegaG, lkj_local
    )
    SS[[i]] <- as.data.frame(SS[[i]])
    SS[[i]]$studyno <- i
    SS[[i]]$obsno <- 1:nrow(SS[[i]])
  }
  do.call("rbind", SS)
}

