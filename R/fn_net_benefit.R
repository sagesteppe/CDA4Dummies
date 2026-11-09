#' Net benefit of thresholded predictions at a target prevalence
#'
#' Decision-curve net benefit of treating (surveying) every unit whose predicted
#' probability is at least the threshold \code{pt}, case-weighted to a target
#' landscape prevalence with \code{\link{reweight}}:
#' \deqn{NB(p_t) = TP_w - FP_w \cdot \frac{p_t}{1 - p_t},}
#' where \eqn{TP_w} and \eqn{FP_w} are the prevalence-weighted true- and
#' false-positive fractions. This is the calibration-\emph{dependent} (per-point)
#' route: it thresholds the probabilities as given.
#'
#' @param truth Outcomes coded 1 (presence) / 0 (absence).
#' @param prob Predicted probabilities, same length as \code{truth}.
#' @param pi Target landscape prevalence in (0, 1).
#' @param pt Numeric vector of threshold probabilities.
#'
#' @return Numeric vector of net benefit, one per element of \code{pt}.
#' @seealso \code{\link{nb_binormal}} (best-achievable from an ROC),
#'   \code{\link{nb_treat_all}} (the treat-all reference).
#' @export
nb_exact <- function(truth, prob, pi, pt) {
  w <- reweight(truth, pi)
  n <- sum(w)
  vapply(pt, function(t) {
    pp <- prob >= t
    sum(w[pp & truth == 1]) / n - sum(w[pp & truth == 0]) / n * (t / (1 - t))
  }, numeric(1))
}

#' Best-achievable net benefit from an equal-variance binormal ROC
#'
#' Net benefit read off an equal-variance binormal ROC (separation \code{mu}) at
#' the operating point whose slope equals the cost ratio
#' \code{((1 - pi) / pi) * (pt / (1 - pt))}. Because it always uses the optimal
#' cutoff, this is the calibration-\emph{independent}, best-achievable net benefit
#' for that ROC -- the summary-metric counterpart of \code{\link{nb_exact}}.
#'
#' @param mu Binormal separation (>= 0).
#' @param pi Target landscape prevalence in (0, 1).
#' @param pt Numeric vector of threshold probabilities.
#'
#' @return Numeric vector of net benefit, one per element of \code{pt}.
#' @export
nb_binormal <- function(mu, pi, pt) {
  vapply(pt, function(t) {
    m  <- ((1 - pi) / pi) * (t / (1 - t))   # target ROC slope for this threshold
    cc <- (log(m) + mu^2 / 2) / mu          # score cutoff giving that slope
    pi * pnorm(mu - cc) - (1 - pi) * pnorm(-cc) * (t / (1 - t))
  }, numeric(1))
}

#' Net benefit of the treat-all (survey-everything) reference strategy
#'
#' The decision-curve reference that acts on every unit:
#' \code{pi - (1 - pi) * (pt / (1 - pt))}. Treat-none is identically 0.
#'
#' @param pi Prevalence in (0, 1).
#' @param pt Numeric vector of threshold probabilities.
#'
#' @return Numeric vector of net benefit, one per element of \code{pt}.
#' @export
nb_treat_all <- function(pi, pt) {
  pi - (1 - pi) * (pt / (1 - pt))
}
