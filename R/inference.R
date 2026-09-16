## this file contains wrappers for stan models

##' MCMC sampling for data from single study
##'
##' This generates a Stan sample for MVN data from a single study
##'
##' @title mvn_infer
##' @param Y multivariate response data - each line sign/symptom values for a patient
##' @param X multivariable explanatory data - each line a set of predictors for a patient
##' @param Z optional new covariate data to predict/simulate Y for (NP columns, same as X;
##'   any number of rows). If supplied, the returned stanfit also carries a generated
##'   quantity `Ynew` (NewObs x NV) with one posterior-predictive draw of Y per MCMC
##'   iteration for each row of Z - extract with \code{mvn_extract_predictions()}.
##'   Default NULL fits the model with no prediction (NewObs=0).
##' @param beta_prior_sd prior for Betas, default=5
##' @param tau_prior_sd = 2.5,
##' @param prior for tau Cauchy scale in correlations, default=2.5
##' @param lkj_prior_scale LKJ prior scale for correlation, default=2
##' @param iter iterations for MCMC, default = 2e3
##' @param cores number of cores to use
##' @param chains number of chains to use
##' @param ...
##' @return a Stan sample object
##' @author Pete Dodd
##' @import rstan
##' @export
mvn_infer <- function(Y, X, Z = NULL,
                      beta_prior_sd = 5, # prior for Betas, def=5
                      tau_prior_sd = 2.5, # prior for tau,   def=2.5
                      lkj_prior_scale = 2, # prior for cor, def=2
                      iter = 2e3, cores = 4, chains = 4, ...) {
  if (!is.null(Z) && ncol(Z) != ncol(X)) {
    stop(
      "Z must have the same number of columns as X (",
      ncol(X), "), got ", ncol(Z)
    )
  }
  ## prepare data
  sdata <- list(
    Nobs = nrow(Y), # number of observations
    NP = ncol(X), # number of variables
    NV = ncol(Y), # number of variates
    X = X, # covariate data
    Y = Y, # outcomes
    beta_prior_sd = beta_prior_sd, # prior for Betas, def=5
    tau_prior_sd = tau_prior_sd, # prior for tau,   def=2.5
    lkj_prior_scale = lkj_prior_scale, # prior for cor, def=2
    NewObs = if (is.null(Z)) 0 else nrow(Z), # number of new obs to predict for
    Z = if (is.null(Z)) {
      matrix(numeric(0), 0, ncol(X))
    } else {
      Z
    } # covariate data for prediction
  )

  ## sample
  rstan::sampling(stanmodels$mvn_infer1e,
    data = sdata,
    chains = chains,
    cores = cores,
    iter = iter, ...
  )
}



##' Extract posterior-predictive draws for new covariate data from a mvn_infer() fit
##'
##' Pulls the \code{Ynew} generated quantity (populated only if \code{mvn_infer()} was
##' called with a non-NULL \code{Z}) out of a stanfit and reshapes it into a proper
##' 3-D array indexed [draw, new observation, variate], rather than leaving the caller
##' to reshape rstan's flattened, row-major summary output by hand (see the row-major
##' reshape note in \code{mvn_extract_hyperparams()}/the package TODO list).
##'
##' @title mvn_extract_predictions
##' @param fit a stanfit object returned by \code{mvn_infer()}, called with non-NULL Z
##' @return a 3-D array [ndraws, NewObs, NV] of posterior-predictive draws of Ynew
##' @author Pete Dodd
##' @export
mvn_extract_predictions <- function(fit) {
  ## NOTE: calling rstan::extract(fit, pars = "Ynew") directly throws
  ## ("no parameter Ynew") rather than returning NULL when Ynew has zero
  ## rows (NewObs=0)
  all_draws <- rstan::extract(fit)
  if (!("Ynew" %in% names(all_draws)) || prod(dim(all_draws$Ynew)[-1]) == 0) {
    stop("fit has no (non-empty) Ynew: was mvn_infer() called with a non-NULL Z?")
  }
  all_draws$Ynew # rstan::extract() already returns this correctly shaped as [draw, NewObs, NV]
}




##' MCMC sampling for data from multiple studies
##'
##' This generates a Stan sample for MVN data from a multiple studies
##'
##' @title mvn_infer_mlm
##' @param Y multivariate response data - each line sign/symptom values for a patient
##' @param X multivariable explanatory data - each line a set of predictors for a patient
##' @param study vector of which study each record belongs to
##' @param betaM_prior_sd prior for Betas
##' @param betaS_prior_sd prior for Betas
##' @param tauM_prior_sd prior for tau
##' @param tauS_prior_sd prior for tau
##' @param lkj_local_prior_scale prior for local cor
##' @param lkj_global_prior_scale prior for global cor
##' @param rhoA beta parameter for rho
##' @param rhoB beta parameter for rho
##' @param iter iterations for MCMC, default
##' @param cores number of cores to use
##' @param chains number of chains to use
##' @param ...
##' @return a Stan sample object
##' @author Pete Dodd
##' @export
##' @import rstan
mvn_infer_mlm <- function(Y, X, study,
                          betaM_prior_sd = 1, # prior for Betas
                          betaS_prior_sd = 0.5, # prior for Betas
                          tauM_prior_sd = 1, # prior for tau
                          tauS_prior_sd = 0.5, # prior for tau
                          lkj_local_prior_scale = 3, # prior for local cor
                          lkj_global_prior_scale = 2, # prior for global cor
                          rhoA = 2, # beta parameter for rho
                          rhoB = 2, # beta parameter for rho
                          iter = 2e3, cores = 4, chains = 4, ...) {
  ## prepare data
  shdata <- list(
    Nrecords = nrow(Y), # number of records/patients
    Nstudies = length(unique(study)), # number of distinct studies
    study = study, # which study does each record correspond to? [Nrecords]
    NP = ncol(X), # number of variables
    NV = ncol(Y), # number of variates
    X = X, # covariate data [Nrecords,NP]
    Y = Y, # outcomes [Nrecords,NV]
    betaM_prior_sd = betaM_prior_sd, # prior for Betas
    betaS_prior_sd = betaS_prior_sd, # prior for Betas
    tauM_prior_sd = tauM_prior_sd, # prior for tau
    tauS_prior_sd = tauS_prior_sd, # prior for tau
    lkj_local_prior_scale = lkj_local_prior_scale, # prior for local cor
    lkj_global_prior_scale = lkj_global_prior_scale, # prior for global cor
    rhoA = rhoA, # beta parameter for rho
    rhoB = rhoB # beta parameter for rho
  )

  ## sample
  rstan::sampling(stanmodels$mvn_inferH,
    data = shdata,
    chains = chains,
    cores = cores,
    iter = iter, ...
  )
}


##' MCMC sampling for data from multiple studies, with horseshoe sparsity on
##' the global correlation structure
##'
##' As \code{mvn_infer_mlm()}, but the global correlation
##' \code{Omega_global} gets a horseshoe shrinkage prior on its off-diagonal
##' entries (promoting sparsity as the number of variates grows) instead of
##' a plain LKJ prior. Each study's local correlation \code{Omega_local[i]}
##' uses a plain LKJ prior.
##'
##' @title mvn_infer_mlm_sparse
##' @param Y multivariate response data - each line sign/symptom values for a patient
##' @param X multivariable explanatory data - each line a set of predictors for a patient
##' @param study vector of which study each record belongs to
##' @param betaM_prior_sd prior for Betas
##' @param betaS_prior_sd prior for Betas
##' @param tauM_prior_sd prior for tau
##' @param tauS_prior_sd prior for tau
##' @param lkj_local_prior_scale prior for local (per-study) correlation
##' @param rhoA beta parameter for rho
##' @param rhoB beta parameter for rho
##' @param iter iterations for MCMC, default
##' @param cores number of cores to use
##' @param chains number of chains to use
##' @param ...
##' @return a Stan sample object
##' @author Pete Dodd
##' @export
##' @import rstan
mvn_infer_mlm_sparse <- function(Y, X, study,
                                 betaM_prior_sd = 1, # prior for Betas
                                 betaS_prior_sd = 0.5, # prior for Betas
                                 tauM_prior_sd = 1, # prior for tau
                                 tauS_prior_sd = 0.5, # prior for tau
                                 lkj_local_prior_scale = 3, # prior for local cor
                                 rhoA = 2, # beta parameter for rho
                                 rhoB = 2, # beta parameter for rho
                                 iter = 2e3, cores = 4, chains = 4, ...) {
  ## prepare data
  shdata <- list(
    Nrecords = nrow(Y), # number of records/patients
    Nstudies = length(unique(study)), # number of distinct studies
    study = study, # which study does each record correspond to? [Nrecords]
    NP = ncol(X), # number of variables
    NV = ncol(Y), # number of variates
    X = X, # covariate data [Nrecords,NP]
    Y = Y, # outcomes [Nrecords,NV]
    betaM_prior_sd = betaM_prior_sd, # prior for Betas
    betaS_prior_sd = betaS_prior_sd, # prior for Betas
    tauM_prior_sd = tauM_prior_sd, # prior for tau
    tauS_prior_sd = tauS_prior_sd, # prior for tau
    lkj_local_prior_scale = lkj_local_prior_scale, # prior for local cor
    rhoA = rhoA, # beta parameter for rho
    rhoB = rhoB # beta parameter for rho
  )

  ## sample
  rstan::sampling(stanmodels$mvn_inferH_sparse,
    data = shdata,
    chains = chains,
    cores = cores,
    iter = iter, ...
  )
}



