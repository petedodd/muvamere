## this file contains MCMC convergence diagnostic helpers for stanfits
## returned by mvn_infer()/mvn_infer_mlm()


##' Summarize MCMC convergence diagnostics for a muvamere stanfit
##'
##' Wraps rstan's standard convergence diagnostics (Rhat, effective sample
##' size, divergent transitions, maximum treedepth hits, E-BFMI) into a
##' single check with an overall pass/fail flag.
##'
##' NOTE on Rhat/ESS: some parameters (e.g. the structural zeros in the
##' upper triangle of a \code{cholesky_factor_corr} such as \code{L_Omega})
##' are constant across every draw by construction, which makes their Rhat
##' and effective sample size come back as \code{NaN} from rstan. This is
##' normal and not a sign of a convergence problem, so those entries are
##' excluded (\code{na.rm = TRUE}) rather than counted as failures.
##'
##' @title mvn_diagnose
##' @param fit a stanfit object, e.g. from
##'   \code{mvn_infer()}/\code{mvn_infer_mlm()}
##' @param rhat_threshold flag as a problem if the max (non-NA) Rhat exceeds
##'   this, default 1.01
##' @param ess_ratio_threshold flag as a problem if the min (non-NA) effective
##'   sample size, as a fraction of total post-warmup draws, falls below this,
##'   default 0.1
##' @param bfmi_threshold flag as a problem if any chain's E-BFMI falls below
##'   this, default 0.2
##' @return an object of class \code{muvamere_diagnostics} (a list) with
##'   elements \code{rhat_max}, \code{ess_min}, \code{ess_min_ratio},
##'   \code{n_divergent}, \code{n_max_treedepth}, \code{bfmi_min},
##'   \code{n_draws}, and \code{ok} (TRUE only if every check passes)
##' @author Pete Dodd
##' @export
mvn_diagnose <- function(fit, rhat_threshold = 1.01, ess_ratio_threshold = 0.1,
                         bfmi_threshold = 0.2) {
  s <- rstan::summary(fit)$summary
  n_draws <- (fit@sim$iter - fit@sim$warmup) * fit@sim$chains

  rhat_max <- suppressWarnings(max(s[, "Rhat"], na.rm = TRUE))
  ess_min <- suppressWarnings(min(s[, "n_eff"], na.rm = TRUE))
  ess_min_ratio <- ess_min / n_draws

  n_divergent <- rstan::get_num_divergent(fit)
  n_max_treedepth <- rstan::get_num_max_treedepth(fit)
  bfmi <- rstan::get_bfmi(fit)
  bfmi_min <- if (length(bfmi)) min(bfmi) else NA_real_

  ok <- is.finite(rhat_max) && rhat_max <= rhat_threshold &&
    is.finite(ess_min_ratio) && ess_min_ratio >= ess_ratio_threshold &&
    n_divergent == 0 &&
    n_max_treedepth == 0 &&
    is.finite(bfmi_min) && bfmi_min >= bfmi_threshold

  structure(
    list(
      rhat_max = rhat_max,
      ess_min = ess_min,
      ess_min_ratio = ess_min_ratio,
      n_divergent = n_divergent,
      n_max_treedepth = n_max_treedepth,
      bfmi_min = bfmi_min,
      n_draws = n_draws,
      ok = ok
    ),
    class = "muvamere_diagnostics"
  )
}

##' @export
print.muvamere_diagnostics <- function(x, ...) {
  cat(if (x$ok) "OK: " else "PROBLEMS FOUND: ", "MCMC diagnostics\n", sep = "")
  cat(sprintf(
    "  Rhat max:        %.4f\n  ESS min:         %.1f (%.1f%% of %d draws)\n  Divergences:     %d\n  Max treedepth:   %d\n  E-BFMI min:      %.3f\n",
    x$rhat_max,
    x$ess_min,
    100 * x$ess_min_ratio, x$n_draws,
    x$n_divergent,
    x$n_max_treedepth,
    x$bfmi_min
  ))
  invisible(x)
}

