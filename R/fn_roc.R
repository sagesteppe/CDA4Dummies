#' ROC-AUC via the rank (Mann-Whitney) identity
#'
#' Area under the ROC curve computed from the rank-sum of the positive scores,
#' equivalent to the probability that a random presence outscores a random
#' absence. Prevalence-invariant, so no case weighting is needed.
#'
#' @param truth Outcomes coded 1 (presence) / 0 (absence).
#' @param prob Predicted probabilities or scores, same length as \code{truth}.
#'
#' @return Scalar ROC-AUC, or \code{NA_real_} if a class is absent.
#' @export
roc_auc <- function(truth, prob) {
  if (length(unique(truth)) < 2) return(NA_real_)
  r  <- rank(prob)
  n1 <- sum(truth == 1)
  n0 <- sum(truth == 0)
  (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
