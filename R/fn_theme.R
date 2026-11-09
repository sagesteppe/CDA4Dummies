#' Shared plotting palette, theme, and scales for the decision-curve figures
#'
#' House style for every figure in the manuscript and appendices, so the look is
#' set in one place. The colours live in \code{R/palette.csv} (role, hex); edit
#' that file to recolour the figures. \code{cda_cols} is that table as a named
#' vector, \code{theme_cda()} is \code{theme_bw()} with the project's base size and
#' legend placement, \code{scale_linetype_reference()} is the dashed-treat-all /
#' solid-treat-none linetype scale, and \code{scale_colour_ratio()} is the natural-
#' prevalence gradient used on the prevalence sweeps.

#' Named colour palette, loaded from R/palette.csv
#'
#' A named character vector mapping semantic roles to hex colours. Roles:
#' \code{treat_none}, \code{treat_all}, \code{naive}, \code{aoa}, \code{summary}
#' (decision-curve series); \code{sweep_point}, \code{sweep_trend} (gap-vs-ratio
#' scatter); \code{ratio_high}, \code{ratio_low} (prevalence-sweep gradient ends).
#' @export
cda_cols <- local({
  path <- file.path(if (exists("HERE", inherits = TRUE)) get("HERE", inherits = TRUE) else "..",
                    "R", "palette.csv")
  tab <- utils::read.csv(path, stringsAsFactors = FALSE)
  stats::setNames(trimws(tab$hex), trimws(tab$role))
})

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
                   plot.caption = ggplot2::element_text(colour = "grey45"))
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

#' Natural-prevalence colour gradient for the prevalence sweeps
#'
#' A continuous (log) gradient running from high P:A ratio / low natural prevalence
#' (\code{ratio_high}) to low ratio / high natural prevalence (\code{ratio_low}),
#' mapped onto a \code{nat_prev} aesthetic.
#'
#' @return A \code{\link[ggplot2]{scale_colour_gradient}} to add to a plot.
#' @export
scale_colour_ratio <- function() {
  ggplot2::scale_colour_gradient(
    name = "natural\nprevalence", trans = "log",
    breaks = c(0.02, 0.05, 0.1, 0.25, 0.5),
    low = cda_cols[["ratio_high"]], high = cda_cols[["ratio_low"]])
}
