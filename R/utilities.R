## this is for functions not relying on Stan



##' MVN sample for a single study
##'
##' MVN sample for a single study using explanatory variables X, regression coefficient matrix Beta, and covariance matrix Sigma
##'
##' @title mvn_simulate
##' @param X multivariable explanatory data - each line a set of predictors for a patient
##' @param Beta matrix of regression coefficients
##' @param Sigma covariance matrix for response
##' @return a matrix of responses Y, with nrow(Y)==nrow(X) and ncol(Y)==ncol(Sigma)
##' @author Pete Dodd
##' @import MASS
##' @export
mvn_simulate <- function(X, Beta, Sigma) {
  MUMatrix <- X %*% Beta
  t(apply(MUMatrix, 1, function(x) MASS::mvrnorm(1, mu = x, Sigma = Sigma)))
}


##' Simulate responses for single study
##'
##' Simulate responses for single study given covariate data.
##' 
##' @title sample_study
##' @param X covariate data
##' @param betag matrix global regression parameters mean
##' @param sigb matrix global regression parameters SD
##' @param rho interpolant between local and global correlations
##' @param taug global tau mean(s)
##' @param sigt global tau SD(s)
##' @param OmegaG global correlation parameters
##' @param lkj_local local correlation LKJ prior parameter
##' @return matrix of responses
##' @author Pete Dodd
##' @export
##' @import trialr
sample_study <- function(X, betag, sigb,
                         rho, taug, sigt,
                         OmegaG, lkj_local) {
  ## dimensions
  NV <- ncol(betag)
  Nobs <- nrow(X)
  NP <- ncol(X) # =nrow(betag)
  ## means
  betas <- matrix(rnorm(prod(dim(betag)), mean = c(betag), sd = c(sigb)),
    nrow = nrow(betag), ncol = ncol(betag)
  )
  ## correlations
  taus <- abs(rnorm(nrow(OmegaG), mean = taug, sd = sigt))
  OmegaL <- trialr::rlkjcorr(1, K = nrow(OmegaG), eta = lkj_local)
  omega <- (rho * OmegaG + (1 - rho) * OmegaL)
  Sig <- diag(taus) %*% omega %*% diag(taus)
  ## samples
  mvn_simulate(X, betas, Sig)
}

##' Simulate responses for a number of studies
##'
##' Given a list of covariate data inputs, a corresponding set of responses is generated.
##'
##' @title simulate_studies
##' @param Xlist list of covariate data matrices 
##' @param betag matrix global regression parameters mean
##' @param sigb matrix global regression parameters SD
##' @param rhoA beta prior parameter for interpolant between local and global correlations
##' @param rhoB beta prior parameter for interpolant between local and global correlations
##' @param taug global tau mean(s)
##' @param sigt global tau SD(s)
##' @param lkj_local local correlation LKJ prior parameter
##' @param lkj_global global correlation LKJ prior parameter
##' @return a data frame with responses for all studies, with studyno and obsno columns appended
##' @author Pete Dodd
##' @import trialr
##' @export
simulate_studies <- function(Xlist,
                             betag, sigb,
                             rhoA, rhoB,
                             taug, sigt,
                             lkj_local, lkj_global) {
  ## sample globals:
  rho <- rbeta(1, rhoA, rhoB)
  OmegaG <- trialr::rlkjcorr(1, nrow(betag), lkj_global)
  ## loop over studies
  SS <- list()
  for (i in 1:length(Xlist)) {
    SS[[i]] <- sample_study(
      Xlist[[i]], betag, sigb,
      rho, taug, sigt,
      OmegaG, lkj_local
    )
    SS[[i]] <- as.data.frame(SS[[i]])
    SS[[i]]$studyno <- i
    SS[[i]]$obsno <- 1:nrow(SS[[i]])
  }
  do.call("rbind", SS)
}
