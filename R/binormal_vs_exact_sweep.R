## ROBUSTNESS SWEEP (exploratory): does the summary-metric (binormal) decision curve agree with
## the exact per-point one, and if the gap varies, does it track anything we can measure?
##
## Per hold-out model (raw, classic excluded) we compute, at the landscape prevalence PI:
##   exact_nb(p_t)   -- net benefit thresholding the RAW probs (prob >= p_t), reweighted to PI.
##                      This is the per-point route; it is calibration-DEPENDENT (uses actual cutoffs).
##   binorm_nb(p_t)  -- net benefit from a binormal ROC fit to the model's reported PR-AUC, read at
##                      the ROC point whose slope matches the cost ratio. This is the summary route;
##                      it is the BEST-achievable NB for that ROC (calibration-INDEPENDENT).
## The gap therefore has two sources: (a) binormal shape approximation, (b) the model's own
## calibration (raw-prob thresholding vs optimal). (b) is why hold-out cox calibration is a natural
## covariate -- NOT a claim, the thing the sweep is meant to test.
##
## Covariates per model: P:A ratio (from name), hold-out prevalence & n, cox intercept & slope.
## Outputs: a per-model table + two exploratory scatterplots + printed correlations. n varies a lot
## across models, so the gap is noisy -- this is a look, not a law.
suppressPackageStartupMessages({library(data.table); library(dplyr); library(ggplot2)})
HERE <- "/home/sagesteppe/Documents/assoRted/CDA4Dummies"
for (f in list.files(file.path(HERE, "R"), "^fn_.*\\.R$", full.names = TRUE)) source(f)
PI   <- 0.025
PT   <- seq(0.005, 0.25, by = 0.005)              # decision-relevant threshold range for the gap
GRID <- 3000                                      # binormal integration grid for this sweep
## helpers: ap_pointwise(), ap_binormal() (fn_average_precision.R);
##   nb_exact() = calibration-dependent exact NB, nb_binormal() = best-achievable NB from the ROC
##   (fn_net_benefit.R); cox_calibration() = hold-out calibration (fn_calibration.R)

d <- fread(file.path(HERE, "data", "holdout_point_predictions.csv"))
d <- d[calibration == "raw" & design == "spatial_mapping"]
d <- d[!(seed %in% c(5001, 5002, 5003))]           # drop the 3-seed replicate experiment; keep 5101-5801 reps
d[, pa_ratio := as.numeric(sub(".*PA1:([0-9.]+)DO.*", "\\1", model))]

res <- d[, {
  if (uniqueN(Presence) < 2 || .N < 30) .SD[0]              # need both classes and a usable n
  else {
    pih <- mean(Presence); pr <- ap_pointwise(Presence, prob, pih)
    mu  <- tryCatch(uniroot(function(m) ap_binormal(m, pih, grid = GRID) - pr, c(1e-3, 8))$root,
                    error = function(e) NA_real_)
    gap <- if (is.na(mu)) rep(NA_real_, length(PT)) else nb_binormal(mu, PI, PT) - nb_exact(Presence, prob, PI, PT)
    cc  <- cox_calibration(Presence, prob)
    data.table(design = design[1], resolution = resolution[1], pa_ratio = pa_ratio[1],
               n = .N, prev_holdout = pih, pr_auc = pr, mu = mu,
               cox_intercept = cc[["cox_intercept"]], cox_slope = cc[["cox_slope"]],
               gap_mean = mean(gap), gap_absmean = mean(abs(gap)))   # gap = binormal - exact
  }
}, by = model]

fwrite(res, file.path(HERE, "data", "binormal_vs_exact_sweep.csv"))

## exploratory views: does the (signed) gap track P:A ratio or hold-out calibration?
plot_gap <- function(xvar, xlab) {
  ggplot(res, aes(.data[[xvar]], gap_mean, size = n)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_point(alpha = 0.5, colour = "#2a78d6") +
    geom_smooth(aes(weight = n), method = "loess", se = TRUE, colour = "#d1495b", linewidth = 0.7, show.legend = FALSE) +
    scale_size_area(max_size = 4) +
    labs(x = xlab, y = "Mean net-benefit gap  (summary - exact)",
         title = sprintf("Summary-vs-exact decision-curve gap at %.1f%% prevalence", 100 * PI),
         subtitle = "Exploratory, one point per hold-out model; point size = hold-out n. Above 0 = summary route optimistic.") +
    theme_bw(base_size = 11)
}
ggsave(file.path(HERE, "figures", "sweep_gap_vs_ratio.png"),
       plot_gap("pa_ratio", "Training P:A ratio (1 : x)") + scale_x_log10(), width = 7, height = 5, dpi = 150)
ggsave(file.path(HERE, "figures", "sweep_gap_vs_cox.png"),
       plot_gap("cox_intercept", "Hold-out cox calibration intercept (0 = calibrated-in-the-large)"),
       width = 7, height = 5, dpi = 150)

cat(sprintf("models summarised: %d\n", nrow(res)))
cat("correlations of signed gap (summary - exact) with:\n")
for (v in c("pa_ratio", "cox_intercept", "cox_slope", "prev_holdout", "pr_auc"))
  cat(sprintf("  %-14s r = % .3f (Spearman % .3f)\n", v,
      cor(res[[v]], res$gap_mean, use = "complete.obs"),
      cor(res[[v]], res$gap_mean, use = "complete.obs", method = "spearman")))
cat("\nwrote data/binormal_vs_exact_sweep.csv, figures/sweep_gap_vs_{ratio,cox}.png\n")
