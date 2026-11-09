#' Exact per-point average precision (PR-AUC) at a target prevalence
#'
#' Average precision computed directly from point-level predictions, case-weighted
#' to a target landscape prevalence with \code{\link{reweight}}. Recall is
#' prevalence-invariant; precision is not, so re-pricing the weights re-prices the
#' PR curve. This is the exact route, used whenever the full \code{(truth, prob)}
#' pairs are available.
#'
#' @param truth Outcomes coded 1 (presence) / 0 (absence).
#' @param prob Predicted probabilities (or any score), same length as \code{truth}.
#' @param pi Target landscape prevalence in (0, 1).
#'
#' @return Scalar average precision, or \code{NA_real_} if a class is absent.
#' @seealso \code{\link{ap_binormal}} for the summary-metric counterpart.
#' @export
ap_pointwise <- function(truth, prob, pi) {
  n1 <- sum(truth == 1)
  n0 <- sum(truth == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  w <- reweight(truth, pi)
  o <- order(prob, decreasing = TRUE)
  truth <- truth[o]
  w <- w[o]
  tpw <- cumsum(w * (truth == 1))
  fpw <- cumsum(w * (truth == 0))
  recall <- tpw / pi                    # total presence weight = pi
  prec   <- tpw / (tpw + fpw)
  sum(c(recall[1], diff(recall)) * prec)
}

#' Average precision of an equal-variance binormal ROC at a prevalence
#'
#' Average precision implied by an equal-variance binormal ROC with separation
#' \code{mu} (so \code{TPR = pnorm(qnorm(FPR) + mu)}), integrated at prevalence
#' \code{pi} by the trapezoid rule. This is the summary-metric route: pair it with
#' \code{\link{ap_scalar_to_prev}} to turn a single reported PR-AUC into a full,
#' prevalence-corrected PR curve shape.
#'
#' @param mu Binormal separation (>= 0). Larger \code{mu} = better discrimination.
#' @param pi Prevalence in (0, 1) at which to price precision.
#' @param grid Number of recall points for the trapezoid integral.
#'
#' @return Scalar average precision.
#' @seealso \code{\link{ap_pointwise}}, \code{\link{ap_scalar_to_prev}}.
#' @export
ap_binormal <- function(mu, pi, grid = 4000) {
  f   <- seq(1e-6, 1 - 1e-6, length.out = grid)  # FPR ascending -> recall ascending
  tpr <- pnorm(qnorm(f) + mu)                     # recall (invariant with prevalence)
  prec <- (pi * tpr) / (pi * tpr + (1 - pi) * f)
  r <- c(0, tpr)
  p <- c(1, prec)
  sum(diff(r) * (head(p, -1) + tail(p, -1)) / 2)  # trapezoid AP
}

#' Correct a scalar PR-AUC from its design prevalence to a target prevalence
#'
#' Given a reported PR-AUC measured at design prevalence \code{pi0}, recover the
#' equal-variance binormal ROC that reproduces it (\code{\link{ap_binormal}}), then
#' re-integrate precision at \code{pi_target}. Use when only a summary PR-AUC and
#' the prevalence it was measured at are available.
#'
#' @param pr_auc Reported PR-AUC (average precision) at the design prevalence.
#' @param pi0 Design/evaluation prevalence the \code{pr_auc} was measured at.
#' @param pi_target Target landscape prevalence to re-price to.
#' @param grid Number of recall points passed to \code{\link{ap_binormal}}.
#'
#' @return List with \code{ap_target} (PR-AUC at \code{pi_target}), the fitted
#'   \code{mu}, and the \code{implied_roc_auc} = \code{pnorm(mu / sqrt(2))}.
#' @export
ap_scalar_to_prev <- function(pr_auc, pi0, pi_target, grid = 20000) {
  mu <- uniroot(function(m) ap_binormal(m, pi0, grid = grid) - pr_auc, c(1e-3, 8))$root
  list(ap_target = ap_binormal(mu, pi_target, grid = grid),
       mu = mu,
       implied_roc_auc = pnorm(mu / sqrt(2)))
}

#' Exact precision-recall curve at a target prevalence
#'
#' Full PR curve from point-level predictions, case-weighted to prevalence
#' \code{pi} with \code{\link{reweight}}. The point-level analogue of
#' \code{\link{ap_pointwise}} (which integrates this curve).
#'
#' @param truth Outcomes coded 1 / 0.
#' @param prob Predicted probabilities/scores.
#' @param pi Target prevalence in (0, 1).
#'
#' @return \code{data.frame} with columns \code{recall} and \code{precision},
#'   ordered by descending score.
#' @export
pr_curve_exact <- function(truth, prob, pi) {
  w <- reweight(truth, pi)
  o <- order(prob, decreasing = TRUE)
  truth <- truth[o]
  w <- w[o]
  tpw <- cumsum(w * (truth == 1))
  fpw <- cumsum(w * (truth == 0))
  data.frame(recall = tpw / pi, precision = tpw / (tpw + fpw))
}

#' Reconstructed precision-recall curve from a binormal ROC at a prevalence
#'
#' PR curve implied by an equal-variance binormal ROC (separation \code{mu}) at
#' prevalence \code{pi}. The point-level analogue of \code{\link{ap_binormal}} and
#' the summary-metric counterpart of \code{\link{pr_curve_exact}}.
#'
#' @param mu Binormal separation (>= 0).
#' @param pi Prevalence in (0, 1).
#' @param grid Number of recall points on the curve.
#'
#' @return \code{data.frame} with columns \code{recall} and \code{precision}.
#' @export
pr_curve_binormal <- function(mu, pi, grid = 400) {
  r <- seq(1e-3, 1 - 1e-3, length.out = grid)
  fpr <- pnorm(qnorm(r) - mu)              # invert TPR = pnorm(qnorm(FPR) + mu)
  o <- pi / (1 - pi)
  data.frame(recall = r, precision = (o * r) / (o * r + fpr))
}
