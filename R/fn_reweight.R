#' Prevalence case weights
#'
#' Case weights that re-price an evaluation sample to a target landscape
#' prevalence. Every presence receives weight \code{pi / n1} and every absence
#' \code{(1 - pi) / n0}, so the total presence weight is \code{pi}, the total
#' absence weight is \code{1 - pi}, and the weights sum to 1. Thresholded
#' summaries (recall, precision, net benefit) computed with these weights behave
#' as if the sample had been drawn at prevalence \code{pi}.
#'
#' @param truth Integer/numeric vector of outcomes coded 1 (presence) / 0 (absence).
#' @param pi Target landscape prevalence in (0, 1).
#'
#' @return Numeric vector of weights, the same length as \code{truth}.
#' @examples
#' w <- reweight(c(1, 0, 0, 1), pi = 0.05)
#' sum(w)            # 1
#' sum(w[c(1, 0, 0, 1) == 1])  # 0.05
#' @export
reweight <- function(truth, pi) {
  n1 <- sum(truth == 1)
  n0 <- sum(truth == 0)
  ifelse(truth == 1, pi / n1, (1 - pi) / n0)
}
