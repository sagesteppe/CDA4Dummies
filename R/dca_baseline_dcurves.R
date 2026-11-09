## Baseline decision curve the standard way: the dcurves package (Vickers et al.), so the plot
## carries the canonical Treat-all / Treat-none references. KISS -- one model, one prevalence.
## The hold-out sample is enriched (~0.52 presence), so we pass the landscape prevalence via
## dca()'s `prevalence=` (its case-control adjustment) rather than reweighting by hand.
suppressPackageStartupMessages({library(dcurves); library(dplyr); library(ggplot2)})
HERE <- "/home/sagesteppe/Documents/assoRted/CDA4Dummies"
PI   <- 0.025
EX   <- "1arc-Iteration1-PA1:30DO:0-Seed5101"   # spatial-mapping, 1:30 P:A (hold-out prev ~3.3%)

h <- read.csv(file.path(HERE, "data", "holdout_point_predictions.csv")) |>
  filter(model == EX, calibration == "raw", design == "spatial_mapping")

PT_MARK <- 0.20                                   # "happy to survey 5 sites to find 1" -> p_t = 1/5
PT_GRID <- seq(0.001, 0.50, by = 0.005)
d <- dca(Presence ~ prob, data = h, prevalence = PI,
         thresholds = PT_GRID,                     # run high so the SDM decay is visible (illustrative)
         label = list(prob = "SDM"))

## ---- SUMMARY-METRIC route: reconstruct the decision curve from the reported PR-AUC alone ----
## Fit an equal-variance binormal ROC (mu) that reproduces the model's PR-AUC at its design
## prevalence, then read net benefit off it: at threshold p_t you operate at the ROC point whose
## slope equals the cost ratio ((1-pi)/pi)*(p_t/(1-p_t)) -- analytic for a binormal.
ap_pointwise <- function(truth, prob, pi) {
  n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  w <- ifelse(truth == 1, pi / n1, (1 - pi) / n0)
  o <- order(prob, decreasing = TRUE); truth <- truth[o]; w <- w[o]
  tpw <- cumsum(w * (truth == 1)); fpw <- cumsum(w * (truth == 0))
  sum(c(tpw[1] / pi, diff(tpw / pi)) * (tpw / (tpw + fpw)))
}
ap_binormal <- function(mu, pi, grid = 4000) {
  f <- seq(1e-6, 1 - 1e-6, length.out = grid); tpr <- pnorm(qnorm(f) + mu)
  prec <- (pi * tpr) / (pi * tpr + (1 - pi) * f)
  sum(diff(c(0, tpr)) * (head(c(1, prec), -1) + tail(c(1, prec), -1)) / 2)
}
pi_hold  <- mean(h$Presence)
pr_hold  <- ap_pointwise(h$Presence, h$prob, pi_hold)      # the number a paper would report
mu_fit   <- uniroot(function(m) ap_binormal(m, pi_hold) - pr_hold, c(1e-3, 8))$root
nb_binormal <- sapply(PT_GRID, function(pt) {
  m <- ((1 - PI) / PI) * (pt / (1 - pt))                    # target ROC slope for this threshold
  cc <- (log(m) + mu_fit^2 / 2) / mu_fit                    # score cutoff giving that slope
  fpr <- pnorm(-cc); tpr <- pnorm(mu_fit - cc)
  PI * tpr - (1 - PI) * fpr * (pt / (1 - pt))
})

## dcurves does the net-benefit math; we render it so Treat-all stays visible where the model
## coincides with it (the raw SDM flags every unit at these thresholds -> SDM == Treat All).
lv <- c("SDM (per-point, exact)", "SDM (summary, binormal)",
        "Survey all (treat-all)", "Survey none (treat-none)")
nb <- as_tibble(d) |>
  transmute(strategy = recode(as.character(label),
                              "Treat All" = "Survey all (treat-all)",
                              "Treat None" = "Survey none (treat-none)",
                              "SDM" = "SDM (per-point, exact)"),
            threshold, net_benefit) |>
  bind_rows(tibble(strategy = "SDM (summary, binormal)",
                   threshold = PT_GRID, net_benefit = nb_binormal)) |>
  mutate(strategy = factor(strategy, levels = lv))
is_ref  <- grepl("treat", nb$strategy)
nb_mark <- nb$net_benefit[nb$strategy == "SDM (per-point, exact)"][
             which.min(abs(nb$threshold[nb$strategy == "SDM (per-point, exact)"] - PT_MARK))]

p <- ggplot(mapping = aes(threshold, net_benefit, colour = strategy)) +
  geom_vline(xintercept = PT_MARK, colour = "grey55", linetype = "dotted") +
  annotate("text", x = PT_MARK, y = PI * 1.08, hjust = -0.04, size = 3, colour = "grey35",
           label = sprintf("p[t]==0.20~': survey 5 to find 1'"), parse = TRUE) +
  annotate("point", x = PT_MARK, y = nb_mark, colour = "#2a78d6", size = 2) +
  geom_line(data = nb[is_ref, ], aes(linetype = strategy), linewidth = 1.1) +
  geom_line(data = nb[!is_ref, ], linewidth = 0.8) +
  scale_colour_manual(values = c("SDM (per-point, exact)" = "#2a78d6",
                                 "SDM (summary, binormal)" = "#e8963a",
                                 "Survey all (treat-all)" = "#d1495b",
                                 "Survey none (treat-none)" = "grey45"), name = NULL) +
  scale_linetype_manual(values = c("Survey all (treat-all)" = "12",
                                   "Survey none (treat-none)" = "solid"), guide = "none") +
  scale_x_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(-0.02, PI * 1.12)) +    # clip the deep treat-all tail, keep its zero-crossing
  labs(x = expression("Threshold probability  " * p[t]), y = "Net benefit",
       title = sprintf("Decision curve for one SDM at %.1f%% landscape prevalence", 100 * PI),
       subtitle = "Per-point (raw predictions) vs the summary route (binormal from reported PR-AUC).") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

ggsave(file.path(HERE, "figures", "dca_baseline_dcurves.png"), p, width = 9, height = 5, dpi = 150)
cat(sprintf("model PR-AUC @ holdout prev %.3f = %.3f ; fitted binormal mu = %.2f\n", pi_hold, pr_hold, mu_fit))
cat("wrote figures/dca_baseline_dcurves.png\n")
cat(sprintf("SDM net benefit at p_t = %.2f (survey 5 to find 1) = %.4f\n", PT_MARK, nb_mark))
{ sdm <- nb[nb$strategy == "SDM (per-point, exact)", ]
  tgt <- c(0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50)
  sdm[sapply(tgt, function(t) which.min(abs(sdm$threshold - t))), ] |>
    transmute(threshold = round(threshold, 3), net_benefit = round(net_benefit, 4)) |>
    as.data.frame() |> print(row.names = FALSE) }
