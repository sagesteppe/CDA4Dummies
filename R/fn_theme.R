#' Shared plotting palette, theme, and scales for the decision-curve figures
#'
#' House style for every figure in the manuscript and appendices, so the look is
#' set in one place. `cda_cols` names the recurring colours by role; `theme_cda()`
#' is `theme_bw()` with the project's base size and legend placement; and
#' `scale_linetype_reference()` is the dashed treat-all / solid treat-none scale
#' used on the decision-curve reference strategies.

#' Named colour palette (semantic roles)
#'
#' \describe{
#'   \item{blue}{the SDM / per-point / in-domain curve (the estimate we trust)}
#'   \item{orange}{the summary-metric / binormal / naive route}
#'   \item{red}{treat-all (survey everything), and the exact-vs-summary contrast}
#'   \item{grey}{treat-none (survey nothing)}
#' }
#' @export
cda_cols <- c(
  blue   = "#2a78d6",
  orange = "#e8963a",
  red    = "#d1495b",
  grey   = "grey45"
)

#' Project ggplot2 theme
#'
#' \code{theme_bw()} at the project base size, with legend placement and a muted
#' caption colour folded in. Legend defaults to \code{"bottom"} (the house style for
#' every figure); pass e.g. \code{legend = "right"} or \code{"none"} to override.
#'
#' @param base_size Base font size passed to \code{\link[ggplot2]{theme_bw}}.
#' @param legend Legend position, e.g. \code{"bottom"}, \code{"right"}, \code{"none"}.
#'
#' @return A ggplot2 theme object to add to a plot.
#' @export
theme_cda <- function(base_size = 11, legend = "bottom") {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = legend,
                   plot.caption = ggplot2::element_text(colour = cda_cols[["grey"]]))
}

#' Linetype scale for the decision-curve reference strategies
#'
#' Dashed treat-all, solid treat-none, with no legend key (the strategies are
#' already distinguished by colour).
#'
#' @return A \code{\link[ggplot2]{scale_linetype_manual}} to add to a plot.
#' @export
scale_linetype_reference <- function() {
  ggplot2::scale_linetype_manual(
    values = c("Survey all (treat-all)" = "12", "Survey none (treat-none)" = "solid"),
    guide = "none")
}
