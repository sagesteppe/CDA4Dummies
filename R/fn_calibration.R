#' Clamp probabilities away from 0 and 1
#'
#' Bound probabilities into \code{[eps, 1 - eps]} so that logit transforms stay
#' finite. Used by the calibration helpers before \code{\link[stats]{qlogis}}.
#'
#' @param p Numeric vector of probabilities.
#' @param eps Small margin; values are clamped to \code{[eps, 1 - eps]}.
#'
#' @return Numeric vector the same length as \code{p}.
#' @export
clamp <- function(p, eps = 1e-6) {
  pmin(pmax(p, eps), 1 - eps)
}

#' Cox calibration intercept and slope on held-out predictions
#'
#' Cox (logistic) calibration of predicted probabilities against observed
#' outcomes. The slope comes from \code{glm(truth ~ logit(prob))}; the intercept
#' ("calibration-in-the-large") from \code{glm(truth ~ 1, offset = logit(prob))}.
#' Separated or non-converging fits (|coef| >= 15) are dropped to \code{NA}.
#'
#' @param truth Outcomes coded 1 / 0.
#' @param prob Predicted probabilities.
#'
#' @return Named numeric vector \code{c(cox_intercept, cox_slope)}; either entry
#'   may be \code{NA_real_}.
#' @export
cox_calibration <- function(truth, prob) {
  lp <- qlogis(clamp(prob))
  fit_slope <- tryCatch(glm(truth ~ lp, family = binomial), error = function(e) NULL)
  fit_intc  <- tryCatch(glm(truth ~ 1, offset = lp, family = binomial), error = function(e) NULL)
  ok <- function(f) !is.null(f) && f$converged && all(abs(coef(f)) < 15)
  c(cox_intercept = if (ok(fit_intc)) unname(coef(fit_intc)[1]) else NA_real_,
    cox_slope     = if (ok(fit_slope)) unname(coef(fit_slope)[2]) else NA_real_)
}

#' Recalibrate predicted probabilities to a target prevalence (logit shift)
#'
#' Shift predictions on the logit scale from their observed prevalence
#' \code{p_obs} to a target prevalence \code{pi}, leaving relative ordering intact:
#' \code{plogis(logit(prob) - logit(p_obs) + logit(pi))}. Inputs are clamped with
#' \code{\link{clamp}} to keep the logits finite.
#'
#' @param prob Predicted probabilities.
#' @param p_obs Prevalence the predictions were produced at.
#' @param pi Target landscape prevalence in (0, 1).
#'
#' @return Numeric vector of recalibrated probabilities.
#' @export
recalibrate_to_prevalence <- function(prob, p_obs, pi) {
  plogis(qlogis(clamp(prob)) - qlogis(clamp(p_obs)) + qlogis(pi))
}
