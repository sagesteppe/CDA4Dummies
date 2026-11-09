## Correcting PR-AUC (average precision) from the prevalence it was measured at to a target
## LANDSCAPE prevalence. Recall (TPR) and FPR are prevalence-invariant; precision is not, so
## PR-AUC must be re-priced. Two routes, one target prevalence:
##
##   (A) PER-POINT reweight  -- exact, when you have point-level (truth, prob). Used on the
##       2026 ground-truth points for the example model.
##   (B) SCALAR -> landscape -- when you only have a reported PR-AUC scalar and the prevalence
##       it was measured at (e.g. the test-holdout number). Recover the binormal ROC that
##       reproduces that AP at its design prevalence, then re-integrate precision at the target.
##
## Example model: 1arc-Iteration1-PA1:1.023DO:0-Seed5101  (see data/example_model/README.md)
suppressPackageStartupMessages(library(tidyverse))
HERE <- "/home/sagesteppe/Documents/assoRted/CDA4Dummies"
for (f in list.files(file.path(HERE, "R"), "^fn_.*\\.R$", full.names = TRUE)) source(f)
PI   <- 0.025                                   # target landscape prevalence for this portion

## (A) exact per-point average precision, case-weighted to prevalence pi:  ap_pointwise()
## (B) scalar PR-AUC at design prevalence pi0 -> AP at target prevalence:  ap_scalar_to_prev()
##     (fits a binormal ROC via ap_binormal()).  All: see R/fn_average_precision.R
GRID <- 20000                                   # fine binormal grid for this a-priori portion

## ---- run on the example model --------------------------------------------------------------
pts <- read.csv(file.path(HERE, "data", "example_model", "groundtruth_points.csv"))
pts <- pts[pts$calibration == "raw", ]
pts$domain <- ifelse(pts$SurfaceArea == "Extrapolation", "Extrapolation",
               ifelse(pts$Distance == "<=270m", "Interp <=270m", "Interp >270m"))
pts$domain <- factor(pts$domain, levels = c("Interp <=270m", "Interp >270m", "Extrapolation"))

realized <- pts |>
  group_by(domain) |>
  summarise(n = n(), sample_prev = mean(Presence),
            pr_auc_sample = ap_pointwise(Presence, prob, mean(Presence)),
            pr_auc_landscape = ap_pointwise(Presence, prob, PI), .groups = "drop") |>
  bind_rows(tibble(domain = "ALL", n = nrow(pts), sample_prev = mean(pts$Presence),
                   pr_auc_sample = ap_pointwise(pts$Presence, pts$prob, mean(pts$Presence)),
                   pr_auc_landscape = ap_pointwise(pts$Presence, pts$prob, PI)))

## a priori: correct the reported TEST-HOLDOUT pr_auc (decoupled from landscape) to PI
pr_test <- 0.836                                  # test_holdout_metrics.csv
pi_test <- 1 / (1 + 1.023)                         # PA1:1.023 -> design prevalence ~0.494
apriori <- ap_scalar_to_prev(pr_test, pi_test, PI, grid = GRID)

write.csv(realized, file.path(HERE, "data", "example_model", "pr_auc_realized_by_domain.csv"),
          row.names = FALSE)

cat(sprintf("target landscape prevalence = %.3f (no-skill PR-AUC baseline = prevalence)\n\n", PI))
cat("REALIZED (per-point ground-truth, exact reweight):\n")
print(as.data.frame(realized |> mutate(across(where(is.numeric), ~round(., 3)))), row.names = FALSE)
cat(sprintf("\nA PRIORI (test-holdout pr_auc %.3f at design prev %.3f -> landscape %.3f):\n",
            pr_test, pi_test, PI))
cat(sprintf("  pr_auc @ %.1f%% prevalence   = %.3f   (lift %.1fx over baseline)\n",
            100 * PI, apriori$ap_target, apriori$ap_target / PI))
cat(sprintf("  implied roc_auc             = %.3f   (reported 0.789 -> pr_auc sits optimistic)\n",
            apriori$implied_roc_auc))
cat("\nwrote data/example_model/pr_auc_realized_by_domain.csv\n")

## ============================================================================================
## SIDE-BY-SIDE: summary-metric (reconstructed) vs per-point (exact) PR curve at landscape pi.
## Both start from the SAME model on the SAME hold-out set and target prevalence PI; they differ
## only in the information used -- a scalar PR-AUC + a binormal ROC-shape assumption, versus the
## full set of held-out predictions. Reconstructed and exact are anchored to the same (PR-AUC,
## prevalence) so the gap is purely the shape assumption. classic design excluded (bias).
## ============================================================================================
suppressPackageStartupMessages(library(data.table))
EX <- "1arc-Iteration1-PA1:1.023DO:0-Seed5101"
h  <- fread(file.path(HERE, "data", "holdout_point_predictions.csv"))
h  <- h[model == EX & calibration == "raw" & design != "classic"]
pi_holdout <- mean(h$Presence)                    # prevalence the reported PR-AUC is measured at
O  <- PI / (1 - PI)

## PR curves exact vs reconstructed: pr_curve_exact(), pr_curve_binormal() -- see fn_average_precision.R
ap_exact_hold <- ap_pointwise(h$Presence, h$prob, pi_holdout)   # == reported scalar
mu_ex <- uniroot(function(m) ap_binormal(m, pi_holdout, grid = GRID) - ap_exact_hold, c(1e-3, 8))$root
ap_exact_land <- ap_pointwise(h$Presence, h$prob, PI)
ap_recon_land <- ap_binormal(mu_ex, PI, grid = GRID)

curve <- rbind(
  data.frame(pr_curve_exact(h$Presence, h$prob, PI),
             route = sprintf("Per-point, exact  (PR-AUC = %.2f)", ap_exact_land)),
  data.frame(pr_curve_binormal(mu_ex, PI),
             route = sprintf("Summary metric, reconstructed  (PR-AUC = %.2f)", ap_recon_land))
)

p2 <- ggplot(curve, aes(recall, precision, colour = route)) +
  geom_hline(yintercept = PI, linetype = "dashed", colour = "grey55") +
  annotate("text", x = 0.02, y = PI, label = sprintf("no-skill = prevalence (%.3f)", PI),
           hjust = 0, vjust = -0.5, size = 3, colour = "grey45") +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("#2a78d6", "#d1495b"), name = NULL) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = "Recall (true-positive rate)", y = sprintf("Precision at landscape prevalence %.1f%%", 100 * PI),
       title = "Same model, same target prevalence: what the summary route misses",
       subtitle = "Reconstruction undershoots the exact curve at low recall — the high-precision region a survey operates in.",
       caption = sprintf("Example model %s; classic design excluded.", EX)) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", plot.caption = element_text(colour = "grey45"))

ggsave(file.path(HERE, "figures", "pr_summary_vs_perpoint.png"), p2, width = 9, height = 5, dpi = 150)

cat(sprintf("\nHOLD-OUT exact vs reconstructed (example model, classic excluded, n=%d, holdout prev=%.3f):\n",
            nrow(h), pi_holdout))
cat(sprintf("  exact PR-AUC @ holdout prev  = %.3f  (reproduces reported scalar)\n", ap_exact_hold))
cat(sprintf("  exact PR-AUC @ %.1f%%          = %.3f  (lift %.1fx)\n", 100*PI, ap_exact_land, ap_exact_land/PI))
cat(sprintf("  reconstructed PR-AUC @ %.1f%%  = %.3f  (lift %.1fx)  -- undershoots by %.0f%%\n",
            100*PI, ap_recon_land, ap_recon_land/PI, 100*(ap_exact_land-ap_recon_land)/ap_exact_land))
cat("\nwrote figures/pr_summary_vs_perpoint.png\n")
